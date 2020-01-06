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

float2 f2intersection_delta(float2 a0, float2 a1, float2 b0, float2 b10){
    float2 d_a = a1 - a0;
    float2 r_a0b0;
    int swap_a;
    float d_a_dot_d_b = dot(d_a,b10);
    if(d_a_dot_d_b<0.0f){
        d_a *= -1.0f;
        d_a_dot_d_b *= -1.0f;
        r_a0b0 = a1 - b0;
        swap_a = 1;
    }else{
        r_a0b0 = a0 - b0;
        swap_a = 0;
    }
    float d_a_dot_d_a = dot(d_a,d_a);
    float d_b_dot_d_b = dot(b10,b10);
    float d_ab_dot_d_a = dot(r_a0b0,d_a);
    float d_ab_dot_d_b = dot(r_a0b0,b10);
    float t_num_p = d_a_dot_d_b*d_ab_dot_d_b;
    float t_num_m = d_b_dot_d_b*d_ab_dot_d_a;
    float det = d_a_dot_d_a*d_b_dot_d_b - d_a_dot_d_b*d_a_dot_d_b;
    float t = (t_num_p-t_num_m)/det;
    if(isnan(t)||isinf(t)){
        t = 0.5f;
    }
    if(t<0.0f)t=0.0f;
    if(t>1.0f)t=1.0f;
    if(swap_a){
        return a1 + d_a*t;
    }else{
        return a0 + d_a*t;
    }
}

#define V0_BIT 0b0001
#define V1_BIT 0b0010
#define V2_BIT 0b0100
#define V3_BIT 0b1000
#define GRID_SIZE 16

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

//void PolygonAddVertex(Polygon *p, float2 v){
//    p->vertices[p->N] = v;
//    p->N++;
//}

#define PolygonAddVertex(polygon_ptr,vertex) \
    (polygon_ptr)->vertices[(polygon_ptr)->N++] = vertex

typedef struct {
    uchar code;
    uchar inside_ends[2];
    uchar inside_edge[2];
    uchar vflag_edge[2];
    float2 v_ends[2];
    float2 v_edge[2];
} PixelEdge;

typedef struct {
    float2 v;
    uchar inside;
} PixelVertex;

typedef struct {
    float2 v0;
    float2 v10;
    float2 N;
} SrcVertex;

typedef struct {
    SrcVertex vertices[4];
} SrcPolygon;

uchar f2BisectSrcPolygon(SrcPolygon *sp, float2 v){
    uchar r=0;
    uchar inside_bit = 1;
    for(int e=0;e<4;e++,inside_bit<<=1){
        float2 vv0 = v - sp->vertices[e].v0;
        if(dot(vv0,sp->vertices[e].N)>-1e-5f){
            r |= inside_bit;
        }
    }
    return r;
}

void PixelEdgeBisectSrcPolygon(PixelEdge *pe, SrcPolygon *sp){
    // test for all outside of any edge
    if((~pe->inside_ends[0])&(~pe->inside_ends[1])&0b1111){
        return;
    }
    // find the intersecting edges
    uchar intersecting = pe->inside_ends[0]^pe->inside_ends[1];
    uchar edge_bit = 1;
    for(int e=0;e<4;e++,edge_bit<<=1){
        if(!(edge_bit&intersecting)) continue;
        switch(pe->code){
        case 0:
            if(edge_bit&pe->inside_ends[0]){
                // v0 is inside
                pe->code = 1;
                pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[e].v0,sp->vertices[e].v10);
                pe->inside_edge[1] = f2BisectSrcPolygon(sp,pe->v_edge[1]);
                pe->vflag_edge[1] = edge_bit;
            }else{
                // v1 is inside
                pe->code = 2;
                pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[e].v0,sp->vertices[e].v10);
                pe->inside_edge[0] = f2BisectSrcPolygon(sp,pe->v_edge[0]);
                pe->vflag_edge[0] = edge_bit;
            }
            break;
        case 1:
            // test for all outside
            if((~pe->inside_ends[0])&(~pe->inside_edge[1])&edge_bit){
                pe->code = 0;
                return;
            }
            // test for an intersection
            if((pe->inside_ends[0]^pe->inside_edge[1])&edge_bit){
                if(pe->inside_ends[0]&edge_bit){
                    // v0 is inside
                    pe->code = 1;
                    pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_edge[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[1] = f2BisectSrcPolygon(sp,pe->v_edge[1]);
                    pe->vflag_edge[1] = edge_bit;
                }else{
                    // v_edge[1] is inside
                    pe->code = 3;
                    pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_edge[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[0] = f2BisectSrcPolygon(sp,pe->v_edge[0]);
                    pe->vflag_edge[0] = edge_bit;
                }
            }
            break;
        case 2:
            // test for all outside
            if((~pe->inside_ends[1])&(~pe->inside_edge[0])&edge_bit){
                pe->code = 0;
                return;
            }
            // test for intersection with this edge
            if((pe->inside_ends[1]^pe->inside_edge[0])&edge_bit){
                if(pe->inside_ends[1]&edge_bit){
                    // v1 is inside
                    pe->code = 2;
                    pe->v_edge[0] = f2intersection_delta(pe->v_edge[0],pe->v_ends[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[0] = f2BisectSrcPolygon(sp,pe->v_edge[0]);
                    pe->vflag_edge[0] = edge_bit;
                }else{
                    // edge->v[0] is inside
                    pe->code = 3;
                    pe->v_edge[1] = f2intersection_delta(pe->v_edge[0],pe->v_ends[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[1] = f2BisectSrcPolygon(sp,pe->v_edge[1]);
                    pe->vflag_edge[1] = edge_bit;
                }
            }
            break;
        case 3:
            // test for all outside
            if((~pe->inside_edge[0])&(~pe->inside_edge[1])&0b1111){
                pe->code = 0;
                return;
            }
            // test for intersection with this edge
            if((pe->inside_edge[0]^pe->inside_edge[1])&edge_bit){
                if(pe->inside_edge[0]&edge_bit){
                    pe->code = 3;
                    pe->v_edge[1] = f2intersection_delta(pe->v_edge[0],pe->v_edge[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[1] = f2BisectSrcPolygon(sp,pe->v_edge[1]);
                    pe->vflag_edge[1] = edge_bit;
                }else{
                    pe->code = 3;
                    pe->v_edge[0] = f2intersection_delta(pe->v_edge[0],pe->v_edge[1],sp->vertices[e].v0,sp->vertices[e].v10);
                    pe->inside_edge[0] = f2BisectSrcPolygon(sp,pe->v_edge[0]);
                    pe->vflag_edge[0] = edge_bit;
                }
            }
            break;
        }
    }
}

void PixelEdgeBorderBisectSrcPolygon(PixelEdge *pe, SrcPolygon *sp){
    // the only case to bisect a border edge is when there is
    // one intersection and the rest are all inside
    uchar intersecting = pe->inside_ends[0]^pe->inside_ends[1];
    uchar all_inside = pe->inside_ends[0]&pe->inside_ends[1];
    switch(intersecting){
    case 1:
        if(all_inside==0b1110){
            if(pe->inside_ends[0]&intersecting){
                // v0 is inside
                pe->code = 1;
                pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[0].v0,sp->vertices[0].v10);
                pe->vflag_edge[1] = 0b0001;
            }else{
                // v1 is inside
                pe->code = 2;
                pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[0].v0,sp->vertices[0].v10);
                pe->vflag_edge[0] = 0b0001;
            }
        }
        return;
    case 2:
        if(all_inside==0b1101){
            if(pe->inside_ends[0]&intersecting){
                // v0 is inside
                pe->code = 1;
                pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[1].v0,sp->vertices[1].v10);
                pe->vflag_edge[1] = 0b0010;
            }else{
                // v1 is inside
                pe->code = 2;
                pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[1].v0,sp->vertices[1].v10);
                pe->vflag_edge[0] = 0b0010;
            }
        }
        return;
    case 4:
        if(all_inside==0b1011){
            if(pe->inside_ends[0]&intersecting){
                // v0 is inside
                pe->code = 1;
                pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[2].v0,sp->vertices[2].v10);
                pe->vflag_edge[1] = 0b0100;
            }else{
                // v1 is inside
                pe->code = 2;
                pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[2].v0,sp->vertices[2].v10);
                pe->vflag_edge[0] = 0b0100;
            }
        }
        return;
    case 8:
        if(all_inside==0b0111){
            if(pe->inside_ends[0]&intersecting){
                // v0 is inside
                pe->code = 1;
                pe->v_edge[1] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[3].v0,sp->vertices[3].v10);
                pe->vflag_edge[1] = 0b1000;
            }else{
                // v1 is inside
                pe->code = 2;
                pe->v_edge[0] = f2intersection_delta(pe->v_ends[0],pe->v_ends[1],sp->vertices[3].v0,sp->vertices[3].v10);
                pe->vflag_edge[0] = 0b1000;
            }
        }
        return;
    default:
        return;
    }
}

void PolygonAddVertexFlags(int flags, SrcPolygon *sp, Polygon *polygon){
    switch(flags){
    case 0b0000:
        return;
    case 0b0001:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b0010:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b0011:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b0100:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b0101:
        // this case is handled by another part of the program
        return;
    case 0b0110:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b0111:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b1000:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    case 0b1001:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b1010:
        // handled by another part of the program
        // shouldn't happen
        return;
    case 0b1011:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b1100:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    case 0b1101:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b1110:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    case 0b1111:
        // should never happen
        return;
    }
}

void PolygonAddSingleVFlag(Polygon *polygon, uchar vflag, SrcPolygon *sp){
    switch(vflag){
    case 0b0001:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b0010:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b0100:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b1000:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    }
}

void PolygonAddMultiVFlag(Polygon *polygon, uchar vflag, SrcPolygon *sp){
    switch(vflag){
    case 0b0000:
        return;
    case 0b0011:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b0110:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b0111:
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        return;
    case 0b1001:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b1011:
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        return;
    case 0b1100:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    case 0b1101:
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        PolygonAddVertex(polygon,sp->vertices[0].v0);
        return;
    case 0b1110:
        PolygonAddVertex(polygon,sp->vertices[1].v0);
        PolygonAddVertex(polygon,sp->vertices[2].v0);
        PolygonAddVertex(polygon,sp->vertices[3].v0);
        return;
    }
}

void PolygonAddEdgeSingleVertexForward(Polygon *polygon, PixelEdge *edge, uchar vflag, SrcPolygon *sp)
{
    switch(edge->code){
    case 0:
        if(edge->inside_ends[0]==0b1111){
            PolygonAddVertex(polygon,edge->v_ends[0]);
        }
        break;
    case 1:
        PolygonAddVertex(polygon,edge->v_ends[0]);
        PolygonAddVertex(polygon,edge->v_edge[1]);
        break;
    case 2:
        vflag &= edge->vflag_edge[0];
        if(vflag){
            PolygonAddSingleVFlag(polygon, vflag, sp);
        }
        PolygonAddVertex(polygon,edge->v_edge[0]);
        break;
    case 3:
        vflag &= edge->vflag_edge[0];
        if(vflag){
            PolygonAddSingleVFlag(polygon, vflag, sp);
        }
        PolygonAddVertex(polygon,edge->v_edge[0]);
        PolygonAddVertex(polygon,edge->v_edge[1]);
        break;
    }
}

void PolygonAddEdgeSingleVertexReverse(Polygon *polygon, PixelEdge *edge, uchar vflag, SrcPolygon *sp)
{
    switch(edge->code){
    case 0:
        if(edge->inside_ends[1]==0b1111){
            PolygonAddVertex(polygon,edge->v_ends[1]);
        }
        break;
    case 1:
        vflag &= edge->vflag_edge[1];
        if(vflag){
            PolygonAddSingleVFlag(polygon, vflag, sp);
        }
        PolygonAddVertex(polygon,edge->v_edge[1]);
        break;
    case 2:
        PolygonAddVertex(polygon,edge->v_ends[1]);
        PolygonAddVertex(polygon,edge->v_edge[0]);
        break;
    case 3:
        vflag &= edge->vflag_edge[1];
        if(vflag){
            PolygonAddSingleVFlag(polygon, vflag, sp);
        }
        PolygonAddVertex(polygon,edge->v_edge[1]);
        PolygonAddVertex(polygon,edge->v_edge[0]);
        break;
    }
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

float2 f2conform_axis(float2 v){
    float2 v_abs = fabs(v);
    if(v_abs.x>=v_abs.y){
        float tan_theta = v_abs.y/v_abs.x;
        if(tan_theta<1e-3){
            v_abs.y = 0;
        }
    }else{
        float tan_theta = v_abs.x/v_abs.y;
        if(tan_theta<1e-3){
            v_abs.x = 0;
        }
    }
    return v_abs * sign(v);
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
            for(int i_srcy=i2_src_max.y;i_srcy>=i2_src_min.y;i_srcy--){
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

kernel void affine_transform_aa_exp(
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
    float2 f2_dsrcx_raw;
    f2_dsrcx_raw.x = M_inv[0];
    f2_dsrcx_raw.y = M_inv[1];
    float2 f2_dsrcy;
    f2_dsrcy.x = -M_inv[3];
    f2_dsrcy.y = -M_inv[4];
    //
    // conform deltas that are close to the axis to lie
    // on the axis;
    float2 f2_dsrcx = f2conform_axis(f2_dsrcx_raw);
    f2_dsrcy = f2conform_axis(f2_dsrcy);

    PixelVertex pixelVertices[GRID_SIZE+1][GRID_SIZE+1];
    PixelEdge xEdges[GRID_SIZE+1][GRID_SIZE];
    PixelEdge yEdges[GRID_SIZE][GRID_SIZE+1];
    uchar pixelVFlags[GRID_SIZE][GRID_SIZE];
    Polygon polygon;
    SrcPolygon srcPolygon;
    //
    // determine the winding order of the source polygon
    // by calculating the determinant of the inverse
    //
    float det = M_inv[0]*M_inv[4] - M_inv[1]*M_inv[3];
    bool ccw;
    if(det>=0.0f){
        ccw = true;
    }else{
        ccw = false;
    }
    //
    // initialize the parts of the source polygon
    // that don't change
    //
    if(ccw){
        srcPolygon.vertices[0].v10 = f2_dsrcy;
        srcPolygon.vertices[1].v10 = f2_dsrcx;
        srcPolygon.vertices[2].v10 = -f2_dsrcy;
        srcPolygon.vertices[3].v10 = -f2_dsrcx;
        for(int i=0;i<4;i++){
            float2 v10 = srcPolygon.vertices[i].v10;
            srcPolygon.vertices[i].N = normalize((float2)(-v10.y,v10.x));
        }
    }else{
        // source vertices are in a clock wise winding
        //
        // first initialize the vertices in the opposite
        // order
        //
        srcPolygon.vertices[3].v0 = (float2)(0.0f,0.0f);
        srcPolygon.vertices[2].v0 = f2_dsrcy;
        srcPolygon.vertices[1].v0 = f2_dsrcy + f2_dsrcx;
        srcPolygon.vertices[0].v0 = f2_dsrcx;
        //
        // then find the deltas and the normals
        //
        for(int i0=0;i0<4;i0++){
            int i1 = i0 + 1;
            if(i1==4)i1=0;
            float2 v10 = srcPolygon.vertices[i1].v0 - srcPolygon.vertices[i0].v0;
            srcPolygon.vertices[i0].v10 = v10;
            srcPolygon.vertices[i0].N = normalize((float2)(-v10.y,v10.x));
        }
    }
    float total_src_area = f2cross(srcPolygon.vertices[0].v10,srcPolygon.vertices[1].v10);

    int i_x_dst;
    for(i_x_dst=0;i_x_dst<width;i_x_dst++,dst++,f2_src00+=f2_dsrcx_raw){
        float2 f2_src01 = f2_src00 + f2_dsrcy;
        float2 f2_src11 = f2_src00 + f2_dsrcy + f2_dsrcx;
        float2 f2_src10 = f2_src00 + f2_dsrcx;
        int2 i2_src00 = convert_int2_plus(f2_src00);
        int2 i2_src01 = convert_int2_plus(f2_src01);
        int2 i2_src11 = convert_int2_plus(f2_src11);
        int2 i2_src10 = convert_int2_plus(f2_src10);
        int2 i2_src_max = max(max(i2_src00,i2_src01),max(i2_src11,i2_src10));
        int2 i2_src_min = min(min(i2_src00,i2_src01),min(i2_src11,i2_src10));
        int Npixelx = i2_src_max.x - i2_src_min.x + 1;
        int Npixely = i2_src_max.y - i2_src_min.y + 1;
        //
        // handle the trivial case of all source vertices inside one pixel
        //
        if(Npixelx==1 && Npixely == 1){
            if(i2_src00.x>=0 && i2_src00.x<width && i2_src00.y<=0 && i2_src00.y>-height){
                *dst = src_image[-i2_src00.y*stride + i2_src00.x];
            }else{
                *dst = bkgnd_pixel;
            }
            continue;
        }
        int2 i2_v0 = (int2)(i2_src_min.x,i2_src_max.y);
        float2 v0 = convert_float2(i2_v0);
        float4 f4_dst_color = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
        float total_area=0.0f;

        // initialize the srcPolygon
        if(ccw){
            srcPolygon.vertices[0].v0 = f2_src00;
            srcPolygon.vertices[1].v0 = f2_src01;
            srcPolygon.vertices[2].v0 = f2_src11;
            srcPolygon.vertices[3].v0 = f2_src10;
        }else{
            srcPolygon.vertices[3].v0 = f2_src00;
            srcPolygon.vertices[2].v0 = f2_src01;
            srcPolygon.vertices[1].v0 = f2_src11;
            srcPolygon.vertices[0].v0 = f2_src10;
        }

        // initialize the pixelVertices
        PixelVertex *pixelVertex = &pixelVertices[0][0];
        float2 v0y = v0;
        for(int y=0;y<Npixely+1;y++){
            float2 v = v0y;
            int x;
            for(x=0;x<Npixelx+1;x++,pixelVertex++){
                pixelVertex->v = v;
                pixelVertex->inside = f2BisectSrcPolygon(&srcPolygon, v);
                v+=(float2)(1.0f,0.0f);
            }
            // move the pointer to the next line
            pixelVertex+=(GRID_SIZE+1) - x;
            v0y+=(float2)(0.0f,-1.0f);
        }

        // initialize the pixelVFlags
        uchar *pixelVFlag = &pixelVFlags[0][0];
        for(int y=0;y<Npixely;y++){
            int x;
            for(x=0;x<Npixelx;x++,pixelVFlag++){
                *pixelVFlag = 0;
            }
            pixelVFlag += GRID_SIZE - x;
        }

        if(ccw){
            pixelVFlags[i2_v0.y-i2_src00.y][i2_src00.x-i2_v0.x] |= V0_BIT;
            pixelVFlags[i2_v0.y-i2_src01.y][i2_src01.x-i2_v0.x] |= V1_BIT;
            pixelVFlags[i2_v0.y-i2_src11.y][i2_src11.x-i2_v0.x] |= V2_BIT;
            pixelVFlags[i2_v0.y-i2_src10.y][i2_src10.x-i2_v0.x] |= V3_BIT;
        }else{
            pixelVFlags[i2_v0.y-i2_src00.y][i2_src00.x-i2_v0.x] |= V3_BIT;
            pixelVFlags[i2_v0.y-i2_src01.y][i2_src01.x-i2_v0.x] |= V2_BIT;
            pixelVFlags[i2_v0.y-i2_src11.y][i2_src11.x-i2_v0.x] |= V1_BIT;
            pixelVFlags[i2_v0.y-i2_src10.y][i2_src10.x-i2_v0.x] |= V0_BIT;
        }

        // bisect and accumulate the src pixels
        PixelVertex *pixel00 = &pixelVertices[0][0];
        PixelVertex *pixel10 = &pixelVertices[0][1];
        PixelVertex *pixel01 = &pixelVertices[1][0];
        PixelVertex *pixel11 = &pixelVertices[1][1];
        PixelEdge *xEdgeTop = &xEdges[0][0];
        PixelEdge *xEdgeBottom = &xEdges[1][0];
        PixelEdge *yEdgeLeft = &yEdges[0][0];
        PixelEdge *yEdgeRight = &yEdges[0][1];
        pixelVFlag = &pixelVFlags[0][0];
        int x_src;
        int y_src;
        for(int y=0,y_src=i2_v0.y;y<Npixely;y++,y_src--){
            int x;
            for(x=0,x_src=i2_v0.x;x<Npixelx;x++,x_src++,
                pixel00++,pixel01++,pixel10++,pixel11++,xEdgeTop++,
                xEdgeBottom++,yEdgeLeft++,yEdgeRight++,pixelVFlag++){
                //
                // bisect the new edges
                //
                // test for the top of the grid
                //
                if(y==0){
                    xEdgeTop->code = 0;
                    xEdgeTop->v_ends[0] = pixel00->v;
                    xEdgeTop->inside_ends[0] = pixel00->inside;
                    xEdgeTop->v_ends[1] = pixel10->v;
                    xEdgeTop->inside_ends[1] = pixel10->inside;
                    PixelEdgeBorderBisectSrcPolygon(xEdgeTop,&srcPolygon);
                }
                //
                // test for the left most edge
                //
                if(x==0){
                    yEdgeLeft->code = 0;
                    yEdgeLeft->v_ends[0] = pixel00->v;
                    yEdgeLeft->inside_ends[0] = pixel00->inside;
                    yEdgeLeft->v_ends[1] = pixel01->v;
                    yEdgeLeft->inside_ends[1] = pixel01->inside;
                    PixelEdgeBorderBisectSrcPolygon(yEdgeLeft,&srcPolygon);
                }
                //
                // bisect the fresh bottom edge
                //
                xEdgeBottom->code = 0;
                xEdgeBottom->v_ends[0] = pixel01->v;
                xEdgeBottom->inside_ends[0] = pixel01->inside;
                xEdgeBottom->v_ends[1] = pixel11->v;
                xEdgeBottom->inside_ends[1] = pixel11->inside;
                if(y<(Npixely-1)){
                    PixelEdgeBisectSrcPolygon(xEdgeBottom,&srcPolygon);
                }else{
                    PixelEdgeBorderBisectSrcPolygon(xEdgeBottom,&srcPolygon);
                }
                //
                // initialize the right edge
                //
                yEdgeRight->code = 0;
                yEdgeRight->v_ends[0] = pixel10->v;
                yEdgeRight->inside_ends[0] = pixel10->inside;
                yEdgeRight->v_ends[1] = pixel11->v;
                yEdgeRight->inside_ends[1] = pixel11->inside;
                if(x<(Npixelx-1)){
                    PixelEdgeBisectSrcPolygon(yEdgeRight,&srcPolygon);
                }else{
                    PixelEdgeBorderBisectSrcPolygon(yEdgeRight,&srcPolygon);
                }
                //
                // create the polygon for this pixel
                //
                polygon.N = 0;
                if(*pixelVFlag==0b0001 || *pixelVFlag==0b0010
                        || *pixelVFlag==0b0100 || *pixelVFlag==0b1000
                        || *pixelVFlag==0b0101 || *pixelVFlag==0b1010){
                    PolygonAddEdgeSingleVertexForward(&polygon,yEdgeLeft,*pixelVFlag,&srcPolygon);
                    PolygonAddEdgeSingleVertexForward(&polygon,xEdgeBottom,*pixelVFlag,&srcPolygon);
                    PolygonAddEdgeSingleVertexReverse(&polygon,yEdgeRight,*pixelVFlag,&srcPolygon);
                    PolygonAddEdgeSingleVertexReverse(&polygon,xEdgeTop,*pixelVFlag,&srcPolygon);
                }else{
                    bool iv_drawn = false;
                    //
                    // draw the left edge
                    //
                    switch(yEdgeLeft->code){
                    case 0:
                        if(yEdgeLeft->inside_ends[0]==0b1111){
                            PolygonAddVertex(&polygon,yEdgeLeft->v_ends[0]);
                        }else{
                            PolygonAddMultiVFlag(&polygon, *pixelVFlag, &srcPolygon);
                            iv_drawn = true;
                        }
                        break;
                    case 1:
                        PolygonAddVertex(&polygon,yEdgeLeft->v_ends[0]);
                        PolygonAddVertex(&polygon,yEdgeLeft->v_edge[1]);
                        break;
                    case 2:
                        PolygonAddVertex(&polygon,yEdgeLeft->v_edge[0]);
                        break;
                    case 3:
                        PolygonAddVertex(&polygon,yEdgeLeft->v_edge[0]);
                        PolygonAddVertex(&polygon,yEdgeLeft->v_edge[1]);
                        break;
                    }
                    //
                    // draw the bottom edge
                    //
                    switch(xEdgeBottom->code){
                    case 0:
                        if(xEdgeBottom->inside_ends[0]==0b1111){
                            PolygonAddVertex(&polygon,xEdgeBottom->v_ends[0]);
                        }else{
                            if(!iv_drawn){
                                PolygonAddMultiVFlag(&polygon, *pixelVFlag, &srcPolygon);
                                iv_drawn = true;
                            }
                        }
                        break;
                    case 1:
                        PolygonAddVertex(&polygon,xEdgeBottom->v_ends[0]);
                        PolygonAddVertex(&polygon,xEdgeBottom->v_edge[1]);
                        break;
                    case 2:
                        PolygonAddVertex(&polygon,xEdgeBottom->v_edge[0]);
                        break;
                    case 3:
                        PolygonAddVertex(&polygon,xEdgeBottom->v_edge[0]);
                        PolygonAddVertex(&polygon,xEdgeBottom->v_edge[1]);
                        break;
                    }
                    //
                    // draw the right edge - it's in the reverse direction to
                    // the drawing order
                    //
                    switch(yEdgeRight->code){
                    case 0:
                        if(yEdgeRight->inside_ends[1]==0b1111){
                            PolygonAddVertex(&polygon,yEdgeRight->v_ends[1]);
                        }else{
                            if(!iv_drawn){
                                PolygonAddMultiVFlag(&polygon, *pixelVFlag, &srcPolygon);
                                iv_drawn = true;
                            }
                        }
                        break;
                    case 1:
                        PolygonAddVertex(&polygon,yEdgeRight->v_edge[1]);
                        break;
                    case 2:
                        PolygonAddVertex(&polygon,yEdgeRight->v_ends[1]);
                        PolygonAddVertex(&polygon,yEdgeRight->v_edge[0]);
                        break;
                    case 3:
                        PolygonAddVertex(&polygon,yEdgeRight->v_edge[1]);
                        PolygonAddVertex(&polygon,yEdgeRight->v_edge[0]);
                        break;
                    }
                    //
                    // draw the top edge - it's in the reverse direction
                    // as well
                    //
                    switch(xEdgeTop->code){
                    case 0:
                        if(xEdgeTop->inside_ends[1]==0b1111){
                            PolygonAddVertex(&polygon,xEdgeTop->v_ends[1]);
                        }else{
                            if(!iv_drawn){
                                PolygonAddMultiVFlag(&polygon, *pixelVFlag, &srcPolygon);
                                iv_drawn = true;
                            }
                        }
                        break;
                    case 1:
                        PolygonAddVertex(&polygon,xEdgeTop->v_edge[1]);
                        break;
                    case 2:
                        PolygonAddVertex(&polygon,xEdgeTop->v_ends[1]);
                        PolygonAddVertex(&polygon,xEdgeTop->v_edge[0]);
                        break;
                    case 3:
                        PolygonAddVertex(&polygon,xEdgeTop->v_edge[1]);
                        PolygonAddVertex(&polygon,xEdgeTop->v_edge[0]);
                        break;
                    }
                }

                uchar4 src_color;
                if(x_src>=0 && x_src<width && y_src<=0 && y_src>-height){
                    src_color = src_image[-y_src*stride + x_src];
                }else{
                    src_color = bkgnd_pixel;
                }
                float4 f4_src_color = convert_float4(src_color);
                f4_src_color /= 255.0f;
                float area = PolygonArea(&polygon);
                total_area += area;
                f4_dst_color += f4_src_color * area;
            }
            //
            // advance the pointers to the next line
            //
            int offset = (GRID_SIZE+1) - x;
            pixel00 += offset;
            pixel01 += offset;
            pixel10 += offset;
            pixel11 += offset;
            yEdgeLeft += offset;
            yEdgeRight += offset;
            offset = GRID_SIZE - x;
            xEdgeTop += offset;
            xEdgeBottom += offset;
            pixelVFlag += offset;
        }

        f4_dst_color /= total_src_area;

        float area_error = (total_area-total_src_area)/total_src_area;
        if(fabs(area_error)>0.001){
            *dst = (uchar4)(0,255,0,255);
        }else{
            *dst = convert_uchar4_sat(f4_dst_color*255.0f);
        }

    }
}

