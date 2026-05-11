#ifndef _REDBIRD_MATLAB_DIFFUSION_H
#define _REDBIRD_MATLAB_DIFFUSION_H

#ifdef  RB_DEBUG
    #define PRINTF(x)        printf x             /**< enable debugging in CPU mode */
#else
    #define PRINTF(x)
#endif

#define GET_1ST_POINTER(x,y,z)  if(strcmp(name,#y)==0) {x->y=(z*)mxGetData(item);PRINTF(("reading rbm.%s\n",#y));}
#define GET_ONE_POINTER(x,y,z)  else GET_1ST_POINTER(x,y,z)

#define GET_1ST_FIELD(x,y)  if(strcmp(name,#y)==0) {double *val=mxGetPr(item);x->y=val[0];PRINTF(("rbm.%s=%g;\n",#y,(float)(x->y)));}
#define GET_ONE_FIELD(x,y)  else GET_1ST_FIELD(x,y)

#define GET_VEC3_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];\
        PRINTF(("rbm.%s=[%g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z)));}
#define GET_VEC34_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];if(mxGetNumberOfElements(item)==4) u->v.w=val[3];\
        PRINTF(("rbm.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w)));}
#define GET_VEC4_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];u->v.w=val[3];\
        PRINTF(("rbm.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w)));}

#define ABS(a)    ((a)<0?-(a):(a))
#define MIN(a,b)  ((a)>(b)?(a):(b))
#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MEXERROR(a)  rb3_throw_exception(999,a,__FILE__,__LINE__)
#define R_C0       3.335640951981520e-12f  //1/C0 in s/mm
#define EPS0_MM    8.854187817e-15         //vacuum permittivity in F/mm (Helmholtz/MWT)

typedef struct Float4 {
    double x, y, z, w;
} float4;

typedef struct Float3 {
    double x, y, z;
} float3;

typedef struct Integer4 {
    int x, y, z, w;
} int4;

typedef struct Integer3 {
    int x, y, z;
} int3;

typedef struct Medium {
    double mua;
    double mus;
    double g;
    double n;
} medium;

typedef struct RedbirdConfig {
    double omega;
    double reff;
    int helmholtz;       /**< 0: diffusion (DOT); 1: Helmholtz (MWT) */
    medium bulkprop;     /**< bulk medium for Bayliss-Turkel BC; reused medium slots: mua=eps_r, mus=sigma, g=mu0, n=n */
    float3 rbcorigin;    /**< origin of the Bayliss-Turkel radiation boundary */
    float4 srcpos;
    float4 srcdir;
} Config;

typedef struct MMC_mesh {
    int nn;      /**< number of nodes */
    int ne;      /**< number of elements */
    int nf;      /**< number of surface triangles */
    int prop;    /**< number of media */
    int ntype;   /**< number of property indices, i.e. type length */
    int ntot;
    int e0;
    int isreoriented;
    float3* node;/**< node coordinates */
    int3* face;  /**< boundary triangle node indices */
    int4* elem;  /**< tetrahedron node indices */
    int*  type;  /**< element or node-based media index */
    medium* med; /**< optical property of different media */
    double* evol; /**< volume of an element */
    double* area; /**< area of the triangular face */
    float3* facecenter; /**< Nf x 3, face centroids (Helmholtz/MWT) */
    float3* facenormal; /**< Nf x 3, face outward unit normals (Helmholtz/MWT) */
    double* facer;      /**< Nf x 1, distance from face centroid to rbcorigin (Helmholtz/MWT) */
    int* rows;
    int* cols;
    int* idxcount;
    int* idxsum;
} tetmesh;

typedef struct FEMForward {
    tetmesh* mesh;
    double* Ar;
    double* Ai;
    double* Dr;
    double* Di;
    double* deldotdel;
} Forward;

typedef struct FEMJacobian {
    int isnodal;
    int nsd, nsdcol;
    int nn, ne;
    double* Phir;
    double* Phii;
    double* Jmuar;
    double* Jmuai;
    double* Jdr;
    double* Jdi;
    double* sd;
    double* deldotdel;
    double* elemid;
    double* elembary;
    double* relem;
} Jacobian;

void rb_femmatrix_nodal(Config* cfg, tetmesh* mesh, Forward* fem);
void rb_femmatrix_elem (Config* cfg, tetmesh* mesh, Forward* fem);
void rb_femmatrix_helmholtz_nodal(Config* cfg, tetmesh* mesh, Forward* fem);
void rb_femmatrix_helmholtz_elem (Config* cfg, tetmesh* mesh, Forward* fem);
void rb_fem_bc(Config* cfg, tetmesh* mesh, Forward* fem);
void rb_fem_bc_helmholtz(Config* cfg, tetmesh* mesh, Forward* fem);
void rb_deldotdel(Config* cfg, tetmesh* mesh, Forward* fem);
void rb_femjacobian(Config* cfg, tetmesh* mesh, Jacobian* jac);

void mcx_set_field(const mxArray* root, const mxArray* item, int idx, Config* cfg, tetmesh* mesh, Jacobian* jac);
void config_init(Config* cfg);
void mcx_clearcfg(Config* cfg);
void mesh_init(tetmesh* mesh);
void jacobian_init(Jacobian* jac);
void mesh_clear(tetmesh* mesh);
void forward_init(Forward* fem, tetmesh* m);
void rb3_usage();
int sub2ind(tetmesh* mesh, int i, int j);
extern "C" int rb3_throw_exception(const int id, const char* msg, const char* filename, const int linenum);

#endif
