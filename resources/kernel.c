kernel void affine_transformation(
    global const char4 *src_image, // rgba source pointer
    global char4       *dst_image, // rgba destination pointer
    int                 stride,    // scanline length
    int                 width,     // image width
    int                 height,    // image height
    global const float *M_inv)     // affine transform inverse in the
                                   // format of a glm::mat3
{
    char4 bkgnd_pixel = (char4)(0xff,0,0,0xFF);
    int i_y_dst = get_global_id(0);
    global char4 *dst = &dst_image[i_y_dst*stride];
    float f_y_dst = convert_float(-i_y_dst);
    float2 f2_src;
    f2_src.x = M_inv[3]*f_y_dst + M_inv[6];
    f2_src.y = M_inv[4]*f_y_dst + M_inv[7];
    // change in source for change in destination
    //
    float2 f2_dsrc;
    f2_dsrc.x = M_inv[0];
    f2_dsrc.y = M_inv[1];
    int i_x_dst;
    for(i_x_dst=0;i_x_dst<width;i_x_dst++,dst++){
        int2 i2_src = convert_int2(f2_src);
        if(i2_src.x>=0 && i2_src.x<width && i2_src.y<=0 && i2_src.y>-height){
            *dst = src_image[-i2_src.y*stride + i2_src.x];
        }else{
            *dst = bkgnd_pixel;
        }
        f2_src += f2_dsrc;
    }
}

float f2cross(float2 a,float2 b){
    return a.x*b.y - a.y*b.x;
}

float2 f2intersection(float2 a0, float2 a1, float2 b0, float2 b1){
    float2 d_a = a1 - a0;
    float2 d_b = b1 - b0;
    float2 r_a0b0 = a0 - b0;
    float d_a_dot_d_a = dot(d_a,d_a);
    float d_b_dot_d_b = dot(d_b,d_b);
    float d_a_dot_d_b = dot(d_a,d_b);
    float t = (d_a_dot_d_b*dot(r_a0b0,d_b)-d_b_dot_d_b*dot(r_a0b0,d_a))/
            (d_a_dot_d_a*d_b_dot_d_b - d_a_dot_d_b*d_a_dot_d_b);
    if(isnan(t)||isinf(t)){
        t=0.5f;
    }
    if(t<0.0f)t=0.0f;
    if(t>1.0f)t=1.0f;
    return a0 + d_a*t;
}

typedef struct{
    int N;
    float2 vertices[10];
} Polygon;

float PolygonArea(Polygon *p){
    if(p->N<3) return 0.0f;
    int Ntri = p->N - 2;
    float2 v0 = p->vertices[0];
    float area = 0.0f;
    for(int t=0;t<Ntri;t++){
        float2 v1 = p->vertices[t+1];
        float2 v2 = p->vertices[t+2];
        area += f2cross(v1-v0,v2-v1);
    }
    return area/2.0f;
}

void PolygonAddVertex(Polygon *p, float2 v){
    p->vertices[p->N] = v;
    p->N++;
}

void PolygonBisectEdge(Polygon *p, Polygon *pr, float2 v0, float2 v1){
    int inside[16];
    pr->N = 0;
    if(p->N<3){
        return;
    }
    float2 v10 = v1 - v0;
    //
    // the cross product of the z-axis vector
    // and the edge vector v10
    float2 N = (float2)(-v10.y,v10.x);
    int all_inside = 1;
    int all_outside = 1;
    for(int v=0;v<p->N;v++){
        float2 vv0 = p->vertices[v] - v0;
        if(dot(N,vv0)>=0.0f){
            inside[v] = 1;
            all_outside = 0;
        }else{
            inside[v] = 0;
            all_inside = 0;
        }
    }
    if(all_inside){
        for(int v=0;v<p->N;v++){
            pr->vertices[v] = p->vertices[v];
        }
        pr->N = p->N;
        return;
    }
    if(all_outside){
        return;
    }
    for(int i_v0=0;i_v0<p->N;i_v0++){
        int i_v1 = i_v0+1;
        if(i_v1==p->N)i_v1=0;
        if(inside[i_v0]){
            PolygonAddVertex(pr,p->vertices[i_v0]);
            if(!inside[i_v1]){
                PolygonAddVertex(pr,f2intersection(p->vertices[i_v0], p->vertices[i_v1], v0, v1));
            }
        }else{
            if(inside[i_v1]){
                PolygonAddVertex(pr,f2intersection(p->vertices[i_v0], p->vertices[i_v1], v0, v1));
            }
        }
    }
}

Polygon *PolygonBisectPolygon(
        Polygon *p,
        Polygon *pb,
        Polygon *pw1,
        Polygon *pw2){
    Polygon *bisected = NULL;
    Polygon *bisected_next=NULL;
    Polygon *ptemp;
    for(int i_v0=0;i_v0<pb->N;i_v0++){
        int i_v1 = i_v0 + 1;
        if(i_v1==pb->N)i_v1=0;
        if(bisected==NULL){
            bisected = pw1;
            PolygonBisectEdge(p, bisected, pb->vertices[i_v0], pb->vertices[i_v1]);
        }else{
            if(bisected_next==NULL){
                bisected_next=pw2;
            }else{
                ptemp = bisected_next;
                bisected_next = bisected;
                bisected = ptemp;
            }
            PolygonBisectEdge(bisected, bisected_next, pb->vertices[i_v0], pb->vertices[i_v1]);
        }
    }
    return bisected_next;
}

int2 convert_int2_plus(float2 v){
    int2 r = convert_int2(v);
    if(v.x<0.0f) r.x--;
    if(v.y>0.0f) r.y++;
    return r;
}



kernel void affine_transform_aa(
    global const uchar4 *src_image, // rgba source pointer
    global uchar4       *dst_image, // rgba destination pointer
    int                 stride,    // scanline length
    int                 width,     // image width
    int                 height,    // image height
    global const float *M_inv)     // affine transform inverse in the
                                   // format of a glm::mat3
{
    uchar4 bkgnd_pixel = (uchar4)(0,0,0,0);
    int i_y_dst = get_global_id(0);
    global uchar4 *dst = &dst_image[i_y_dst*stride];
    float f_y_dst = convert_float(-i_y_dst);
    float2 f2_src00;
    f2_src00.x = M_inv[3]*f_y_dst + M_inv[6];
    f2_src00.y = M_inv[4]*f_y_dst + M_inv[7];
    // change in source for change in destination
    //
    float2 f2_dsrcx;
    f2_dsrcx.x = M_inv[0];
    f2_dsrcx.y = M_inv[1];
    float2 f2_dsrcy;
    f2_dsrcy.x = -M_inv[3];
    f2_dsrcy.y = -M_inv[4];
    Polygon pixel_polygon;
    Polygon src_polygon;
    Polygon w1_polygon;
    Polygon w2_polygon;
    int i_x_dst;
    for(i_x_dst=0;i_x_dst<width;i_x_dst++,dst++){
        float2 f2_src01 = f2_src00 + f2_dsrcy;
        float2 f2_src11 = f2_src00 + f2_dsrcy + f2_dsrcx;
        float2 f2_src10 = f2_src00 + f2_dsrcx;
        src_polygon.N=4;
        src_polygon.vertices[0] = f2_src00;
        src_polygon.vertices[1] = f2_src01;
        src_polygon.vertices[2] = f2_src11;
        src_polygon.vertices[3] = f2_src10;
        float total_src_area = PolygonArea(&src_polygon);
        int2 i2_src00 = convert_int2_plus(f2_src00);
        int2 i2_src01 = convert_int2_plus(f2_src01);
        int2 i2_src11 = convert_int2_plus(f2_src11);
        int2 i2_src10 = convert_int2_plus(f2_src10);
        int2 i2_src_max = max(max(i2_src00,i2_src01),max(i2_src11,i2_src10));
        int2 i2_src_min = min(min(i2_src00,i2_src01),min(i2_src11,i2_src10));
        float4 f4_dst_color = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
        for(int i_srcx=i2_src_min.x;i_srcx<=i2_src_max.x;i_srcx++){
            for(int i_srcy=i2_src_min.y;i_srcy<=i2_src_max.y;i_srcy++){
                uchar4 src_color;
                if(i_srcx>=0 && i_srcx<width && i_srcy<=0 && i_srcy>-height){
                    src_color = src_image[-i_srcy*stride + i_srcx];
                }else{
                    src_color = bkgnd_pixel;
                }
                float4 f4_src_color = convert_float4(src_color);
                f4_src_color/=255.0f;
                //f4_src_color = (float4)(1.0f,1.0f,1.0f,1.0f);
                pixel_polygon.vertices[0] = convert_float2((int2)(i_srcx, i_srcy));
                pixel_polygon.vertices[1] = convert_float2((int2)(i_srcx, i_srcy-1));
                pixel_polygon.vertices[2] = convert_float2((int2)(i_srcx+1, i_srcy-1));
                pixel_polygon.vertices[3] = convert_float2((int2)(i_srcx+1, i_srcy));
                pixel_polygon.N = 4;
                Polygon *bisected = PolygonBisectPolygon(
                            &pixel_polygon,&src_polygon,
                            &w1_polygon,&w2_polygon);
                float pixel_area = PolygonArea(bisected);
                f4_dst_color += f4_src_color*pixel_area;
            }
        }

        f4_dst_color /= total_src_area;

        *dst = convert_uchar4_sat(f4_dst_color*255.0f);;
        f2_src00 += f2_dsrcx;
    }
}

