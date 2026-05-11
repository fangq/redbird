/*
 * =============================================================
 * FEM Kernel for the diffusion equation (DOT) and the scalar
 * Helmholtz equation (MWT).
 * Author: Qianqian Fang(q.fang < at > neu.edu)
 * Created on: 2007/11/21
 * Version: 0.8.0
 *
 * DOT Model:
 *    [A]*c = b,  [A] = <D del phi_i, del phi_j> + <mua phi_i, phi_j>
 *                  + i*omega/c0*nref <phi_i,phi_j> + Robin BC
 * MWT Model (cfg.helmholtz=1 or cfg.bulkprop given):
 *    [A]*E = b,  [A] = <del phi_i, del phi_j> - <k^2 phi_i, phi_j>
 *                  + 1st-order Bayliss-Turkel radiation BC
 *    where k^2(x) = omega^2*mu0(x)*eps0*eps_r(x) - i*omega*mu0(x)*sigma(x)
 *
 * In MWT mode cfg.prop columns are reinterpreted as
 *    [eps_r, sigma(S/mm), mu0(H/mm), n]
 * and the medium struct slots mua/mus/g/n are reused to carry them.
 * =============================================================
 */

#include <string.h>
#include <exception>
#include <stdio.h>
#include <math.h>
#include <complex>

#include "mex.h"
#include "rbfemmatrix.h"

#define N137_REFF 0.468407140351385

void mexFunction(int nlhs,       mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    Config cfg;
    tetmesh mesh;
    Forward femdata;
    Jacobian jac;
    int isjacobian = 0;
    int nfields, ifield;

    if (nrhs < 1) {
        rb3_usage();
        return;
    }

    /* Get the matrices from the input list */

    if (!mxIsStruct(prhs[0])) {
        mexErrMsgTxt("Input must be a structure.");
    }

    try {
        config_init(&cfg);
        mesh_init(&mesh);
        forward_init(&femdata, &mesh);
        jacobian_init(&jac);

        nfields = mxGetNumberOfFields(prhs[0]);

        for (ifield = 0; ifield < nfields; ifield++) { /* how many input struct fields */
            mxArray* tmp = mxGetFieldByNumber(prhs[0], 0, ifield);

            if (tmp == NULL) {
                continue;
            }

            mcx_set_field(prhs[0], tmp, ifield, &cfg, &mesh, &jac);
        }

        if (mesh.nn <= 0 || mesh.ne <= 0) {
            mexErrMsgTxt("dimensions of the input or output parameters are incorrect.");
        }

        if (nrhs > 1) {
            const mwSize* arraydim;
            arraydim = mxGetDimensions(prhs[1]);

            if (arraydim[1] < 3 || arraydim[0] == 0) {
                MEXERROR("sd must have at least 3 columns");
            }

            jac.sd = mxGetPr(prhs[1]);
            jac.nsd = arraydim[0];
            jac.nsdcol = arraydim[1];
            isjacobian = nlhs;

            if (nrhs <= 2) {
                MEXERROR("sd and phi must be given to calculate Jacobians");
            }

            arraydim = mxGetDimensions(prhs[2]);

            if (arraydim[0] != (const mwSize)mesh.nn) {
                MEXERROR("phi's row number must be the same as the number of node");
            }

            jac.Phir = mxGetPr(prhs[2]); // phi: [mesh.nn x noptodes]

            if (mxIsComplex(prhs[2])) {
                jac.Phii = mxGetPi(prhs[2]);    // phi: [mesh.nn x noptodes]
            }

            if (nrhs > 3) {
                arraydim = mxGetDimensions(prhs[3]);

                if (arraydim[0] != (const mwSize)mesh.ne || arraydim[1] != 10) {
                    MEXERROR("deldotdel must have be a 2D array with size: elem number x 10");
                }

                jac.deldotdel = mxGetPr(prhs[3]); // phi: [mesh.nn x jac.nsd]
            }

            if (nrhs > 4) {
                jac.isnodal = ((int) * (mxGetPr(prhs[4])));
            }

            if (nrhs > 8) {
                jac.isnodal = 1; // build nodal Jmua and Jd on a coarse mesh
                arraydim = mxGetDimensions(prhs[5]);

                if (arraydim[0]*arraydim[1] != (const mwSize)mesh.nn) {
                    MEXERROR("elemid must have the same length has the forward mesh node number");
                }

                jac.elemid = mxGetPr(prhs[5]);

                arraydim = mxGetDimensions(prhs[6]);

                if (arraydim[0] != (const mwSize)mesh.nn && arraydim[1] != 4) {
                    MEXERROR("elembary must have the same lengt has the forward mesh node number and 4 columns");
                }

                jac.elembary = mxGetPr(prhs[6]);
                jac.nn = ((int) * (mxGetPr(prhs[7])));

                if (jac.nn > mesh.nn) {
                    MEXERROR("reconstruction mesh must be coarser than or at most the same as the forward mesh");
                }

                arraydim = mxGetDimensions(prhs[8]);

                if (arraydim[1] != 4) {
                    MEXERROR("relem must have 4 columns");
                }

                jac.ne = arraydim[0];
                jac.relem = mxGetPr(prhs[8]);
            }
        }

        /* In Helmholtz/MWT mode the operator is always complex when omega>0
         * (mass term has +i*omega*mu*sigma) and the Bayliss-Turkel BC has
         * an imaginary part proportional to i*kbg.  Treat any non-zero
         * omega in helmholtz mode the same as the DOT frequency-domain
         * case below.  In the (degenerate) omega=0 helmholtz case the
         * matrix is real and we fall through to the else branch. */
        int wantcomplex = (cfg.omega > 0.f);

        /* Assign a pointer to the output. */
        if (wantcomplex) {
            if (nlhs >= 1) {
                if (isjacobian == 0) {
                    plhs[0] = mxCreateDoubleMatrix(1, mesh.nn, mxCOMPLEX);
                    femdata.Di = mxGetPi(plhs[0]);
                } else {
                    plhs[0] = mxCreateDoubleMatrix(jac.nsd, (jac.nn ? jac.nn : (jac.isnodal ? mesh.nn : mesh.ne)), mxCOMPLEX);
                    jac.Jmuai = mxGetPi(plhs[0]);
                }
            }

            if (nlhs >= 2) {
                if (isjacobian == 0) {
                    plhs[1] = mxCreateDoubleMatrix(1, mesh.ntot, mxCOMPLEX);
                    femdata.Ai = mxGetPi(plhs[1]);
                } else {
                    plhs[1] = mxCreateDoubleMatrix(jac.nsd, (jac.nn ? jac.nn : (jac.isnodal ? mesh.nn : mesh.ne)), mxCOMPLEX);
                    jac.Jdi = mxGetPi(plhs[1]);
                }
            }
        } else {
            if (isjacobian == 0) {
                if (nlhs >= 1) {
                    plhs[0] = mxCreateDoubleMatrix(1, mesh.nn, mxREAL);
                }

                if (nlhs >= 2) {
                    plhs[1] = mxCreateDoubleMatrix(1, mesh.ntot, mxREAL);
                }
            } else {
                if (nlhs >= 1) {
                    plhs[0] = mxCreateDoubleMatrix(jac.nsd, (jac.nn ? jac.nn : (jac.isnodal ? mesh.nn : mesh.ne)), mxREAL);
                }

                if (nlhs >= 2) {
                    plhs[1] = mxCreateDoubleMatrix(jac.nsd, (jac.nn ? jac.nn : (jac.isnodal ? mesh.nn : mesh.ne)), mxREAL);
                }
            }
        }

        if (isjacobian == 0) {
            if (nlhs >= 1) {
                femdata.Dr = mxGetPr(plhs[0]);
            }

            if (nlhs >= 2) {
                femdata.Ar = mxGetPr(plhs[1]);
            }
        } else {
            if (nlhs >= 1) {
                jac.Jmuar = mxGetPr(plhs[0]);
            }

            if (nlhs >= 2) {
                jac.Jdr = mxGetPr(plhs[1]);
            }
        }

        if (nlhs >= 3) {
            plhs[2] = mxCreateDoubleMatrix(10, mesh.ne, mxREAL);
            femdata.deldotdel = mxGetPr(plhs[2]);
        } else {
            femdata.deldotdel = (double*)calloc(mesh.ne, 10 * sizeof(double));
        }

        if (isjacobian == 0) {
            rb_deldotdel(&cfg, &mesh, &femdata);

            if (cfg.helmholtz) {
                if (mesh.ntype == mesh.ne) {
                    rb_femmatrix_helmholtz_elem (&cfg, &mesh, &femdata);
                } else if (mesh.ntype == mesh.nn) {
                    rb_femmatrix_helmholtz_nodal(&cfg, &mesh, &femdata);
                } else {
                    MEXERROR("cfg.seg must have a length that matches the length of either cfg.node or cfg.elem");
                }

                rb_fem_bc_helmholtz(&cfg, &mesh, &femdata);
            } else {
                if (mesh.ntype == mesh.ne) {
                    rb_femmatrix_elem (&cfg, &mesh, &femdata);
                } else if (mesh.ntype == mesh.nn) {
                    rb_femmatrix_nodal(&cfg, &mesh, &femdata);
                } else {
                    MEXERROR("cfg.seg must have a length that matches the length of either cfg.node or cfg.elem");
                }

                rb_fem_bc(&cfg, &mesh, &femdata);
            }
        } else {
            rb_femjacobian(&cfg, &mesh, &jac);
        }

        if (nlhs < 3 && femdata.deldotdel) {
            delete (femdata.deldotdel);
        }

        /** \subsection sclean End the simulation */
        mesh_clear(&mesh);

    } catch (const char* err) {
        mexPrintf("Error: %s\n", err);
    } catch (const std::exception& err) {
        mexPrintf("C++ Error: %s\n", err.what());
    } catch (...) {
        mexPrintf("Unknown Exception");
    }
}

void mcx_set_field(const mxArray* root, const mxArray* item, int idx, Config* cfg, tetmesh* mesh, Jacobian* jac) {
    const char* name = mxGetFieldNameByNumber(root, idx);
    const mwSize* arraydim;
    int i, j;

    GET_1ST_FIELD(cfg, omega)
    GET_VEC34_FIELD(cfg, srcpos)
    GET_VEC34_FIELD(cfg, srcdir)
    GET_VEC3_FIELD(cfg, rbcorigin)
    GET_ONE_FIELD(mesh, e0)
    GET_ONE_FIELD(mesh, isreoriented)
    GET_ONE_FIELD(cfg, reff)
    GET_ONE_FIELD(cfg, helmholtz)
    else if (strcmp(name, "node") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 3) {
            MEXERROR("the 'node' field must have 3 columns (x,y,z)");
        }

        double* val = mxGetPr(item);
        mesh->nn = arraydim[0];

        if (mesh->node) {
            free(mesh->node);
        }

        mesh->node = (float3*)calloc(sizeof(float3), mesh->nn);

        for (j = 0; j < 3; j++)
            for (i = 0; i < mesh->nn; i++) {
                ((double*)(&mesh->node[i]))[j] = val[j * mesh->nn + i];
            }

        PRINTF(("rbm.nn=%d;\n", mesh->nn));
    } else if (strcmp(name, "elem") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 4) {
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        }

        double* val = mxGetPr(item);
        mesh->ne = arraydim[0];

        if (mesh->elem) {
            free(mesh->elem);
        }

        mesh->elem = (int4*)calloc(sizeof(int4), mesh->ne);

        for (j = 0; j < 4; j++)
            for (i = 0; i < mesh->ne; i++) {
                ((int*)(&mesh->elem[i]))[j] = (int)val[j * mesh->ne + i] - 1;
            }

        PRINTF(("rbm.ne=%d;\n", mesh->ne));
    } else if (strcmp(name, "face") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 3) {
            MEXERROR("the 'face' field must have 3 columns (e1,e2,e3)");
        }

        double* val = mxGetPr(item);
        mesh->nf = arraydim[0];

        if (mesh->face) {
            free(mesh->face);
        }

        mesh->face = (int3*)calloc(sizeof(int3), mesh->nf);

        for (j = 0; j < 3; j++)
            for (i = 0; i < mesh->nf; i++) {
                ((int*)(&mesh->face[i]))[j] = (int)val[j * mesh->nf + i] - 1;
            }

        PRINTF(("rbm.nf=%d;\n", mesh->nf));
    } else if (strcmp(name, "seg") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'seg' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ntype = MAX(arraydim[0], arraydim[1]);

        if (mesh->type) {
            free(mesh->type);
        }

        mesh->type = (int*)malloc(sizeof(int ) * mesh->ntype);

        for (i = 0; i < mesh->ntype; i++) {
            mesh->type[i] = val[i];
        }

        PRINTF(("rbm.seg=<%d>;\n", mesh->ntype));
    } else if (strcmp(name, "evol") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'evol' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ne = MAX(arraydim[0], arraydim[1]);

        if (mesh->evol) {
            free(mesh->evol);
        }

        mesh->evol = (double*)malloc(sizeof(double) * mesh->ne);

        for (i = 0; i < mesh->ne; i++) {
            mesh->evol[i] = val[i];
        }

        PRINTF(("rbm.evol=%d;\n", mesh->ne));
    } else if (strcmp(name, "area") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'evol' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->nf = MAX(arraydim[0], arraydim[1]);

        if (mesh->area) {
            free(mesh->area);
        }

        mesh->area = (double*)malloc(sizeof(double) * mesh->nf);

        for (i = 0; i < mesh->nf; i++) {
            mesh->area[i] = val[i];
        }

        PRINTF(("rbm.area=%d;\n", mesh->nf));
    } else if (strcmp(name, "idxcount") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'evol' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->nn = MAX(arraydim[0], arraydim[1]);

        if (mesh->idxcount) {
            free(mesh->idxcount);
        }

        mesh->idxcount = (int*)malloc(sizeof(int) * mesh->nn);

        for (i = 0; i < mesh->nn; i++) {
            mesh->idxcount[i] = val[i];
        }

        PRINTF(("rbm.idxcount=%d;\n", mesh->nn));
    } else if (strcmp(name, "idxsum") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'evol' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->nn = MAX(arraydim[0], arraydim[1]);

        if (mesh->idxsum) {
            free(mesh->idxsum);
        }

        mesh->idxsum = (int*)malloc(sizeof(int) * mesh->nn);

        for (i = 0; i < mesh->nn; i++) {
            mesh->idxsum[i] = val[i];
        }

        PRINTF(("rbm.idxsum=%d;\n", mesh->nn));
    } else if (strcmp(name, "rows") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'rows' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ntot = MAX(arraydim[0], arraydim[1]);

        if (mesh->rows) {
            free(mesh->rows);
        }

        mesh->rows = (int*)malloc(sizeof(int) * mesh->ntot);

        for (i = 0; i < mesh->ntot; i++) {
            mesh->rows[i] = (int)val[i] - 1;
        }

        PRINTF(("rbm.rows=%d;\n", mesh->ntot));
    } else if (strcmp(name, "cols") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'cols' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ntot = MAX(arraydim[0], arraydim[1]);

        if (mesh->cols) {
            free(mesh->cols);
        }

        mesh->cols = (int*)malloc(sizeof(int) * mesh->ntot);

        for (i = 0; i < mesh->ntot; i++) {
            mesh->cols[i] = (int)val[i] - 1;
        }

        PRINTF(("rbm.cols=%d;\n", mesh->ntot));
    } else if (strcmp(name, "prop") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] > 0 && arraydim[1] != 4) {
            MEXERROR("the 'prop' field must have 4 columns (mua,mus,g,n)");
        }

        double* val = mxGetPr(item);
        mesh->prop = arraydim[0] - 1;

        if (mesh->med) {
            free(mesh->med);
        }

        mesh->med = (medium*)calloc(sizeof(medium), mesh->prop + 1);

        for (j = 0; j < 4; j++)
            for (i = 0; i <= mesh->prop; i++) {
                ((double*)(&mesh->med[i]))[j] = val[j * (mesh->prop + 1) + i];
            }

        PRINTF(("rbm.prop=%d;\n", mesh->prop));
    } else if (strcmp(name, "deldotdel") == 0) {
        if (jac->deldotdel == NULL) {
            arraydim = mxGetDimensions(item);

            if (arraydim[1] != 10) {
                MEXERROR("deldotphi must have be a 2D array with size: node number x 10");
            }

            jac->deldotdel = mxGetPr(item); // phi: [mesh.nn x jac.nsd]
            PRINTF(("rbm.deldotdel=<%d,10>;\n", arraydim[0]));
        }
    } else if (strcmp(name, "bulkprop") == 0) {
        /* Bulk medium for the 1st-order Bayliss-Turkel BC (Helmholtz/MWT).
         * Expects a 4-element vector [eps_r, sigma, mu0, n] in matlab order.
         * Implicitly engages helmholtz mode if not already set. */
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) < 3) {
            MEXERROR("the 'bulkprop' field must have at least 3 elements: [eps_r, sigma, mu0, (n)]");
        }

        double* val = mxGetPr(item);
        cfg->bulkprop.mua = val[0];  // eps_r (re-used slot)
        cfg->bulkprop.mus = val[1];  // sigma
        cfg->bulkprop.g   = val[2];  // mu0

        if (MAX(arraydim[0], arraydim[1]) > 3) {
            cfg->bulkprop.n = val[3];
        }

        cfg->helmholtz = 1;
        PRINTF(("rbm.bulkprop=[%g %g %g %g]; helmholtz=1\n", cfg->bulkprop.mua, cfg->bulkprop.mus, cfg->bulkprop.g, cfg->bulkprop.n));
    } else if (strcmp(name, "facecenter") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 3) {
            MEXERROR("the 'facecenter' field must have 3 columns (x,y,z) and Nf rows");
        }

        double* val = mxGetPr(item);
        int n = arraydim[0];

        if (mesh->facecenter) {
            free(mesh->facecenter);
        }

        mesh->facecenter = (float3*)calloc(sizeof(float3), n);

        for (j = 0; j < 3; j++)
            for (i = 0; i < n; i++) {
                ((double*)(&mesh->facecenter[i]))[j] = val[j * n + i];
            }

        PRINTF(("rbm.facecenter=<%d,3>;\n", n));
    } else if (strcmp(name, "facenormal") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 3) {
            MEXERROR("the 'facenormal' field must have 3 columns (nx,ny,nz) and Nf rows");
        }

        double* val = mxGetPr(item);
        int n = arraydim[0];

        if (mesh->facenormal) {
            free(mesh->facenormal);
        }

        mesh->facenormal = (float3*)calloc(sizeof(float3), n);

        for (j = 0; j < 3; j++)
            for (i = 0; i < n; i++) {
                ((double*)(&mesh->facenormal[i]))[j] = val[j * n + i];
            }

        PRINTF(("rbm.facenormal=<%d,3>;\n", n));
    } else if (strcmp(name, "facer") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'facer' field can not be empty");
        }

        double* val = mxGetPr(item);
        int n = MAX(arraydim[0], arraydim[1]);

        if (mesh->facer) {
            free(mesh->facer);
        }

        mesh->facer = (double*)malloc(sizeof(double) * n);

        for (i = 0; i < n; i++) {
            mesh->facer[i] = val[i];
        }

        PRINTF(("rbm.facer=<%d>;\n", n));
    } else if (strcmp(name, "detpos") == 0) {
        // do nothing
    } else if (strcmp(name, "srctype") == 0) {
        // do nothing
    } else if (strcmp(name, "srcparam1") == 0) {
        // do nothing
    } else if (strcmp(name, "srcparam2") == 0) {
        // do nothing
    } else if (strcmp(name, "srcpattern") == 0) {
        // do nothing
    } else if (strcmp(name, "dettype") == 0) {
        // do nothing
    } else if (strcmp(name, "detparam1") == 0) {
        // do nothing
    } else if (strcmp(name, "detparam2") == 0) {
        // do nothing
    } else if (strcmp(name, "detpattern") == 0) {
        // do nothing
    } else if (strcmp(name, "srcpos0") == 0) {
        // do nothing
    } else if (strcmp(name, "detpos0") == 0) {
        // do nothing
    } else if (strcmp(name, "detdir") == 0) {
        // do nothing
    } else if (strcmp(name, "nvol") == 0) {
        // do nothing
    } else if (strcmp(name, "mua") == 0) {
        // do nothing
    } else if (strcmp(name, "musp0") == 0) {
        // do nothing
    } else {
        printf("WARNING: redundant field '%s'\n", name);
    }

    mexEvalString("pause(.0001);");
}

void config_init(Config* cfg) {
    cfg->omega = 0.f;
    cfg->reff = N137_REFF;
    cfg->helmholtz = 0;
    memset(&(cfg->bulkprop), 0, sizeof(medium));
    memset(&(cfg->rbcorigin), 0, sizeof(float3));
    memset(&(cfg->srcpos), 0, sizeof(float4));
    memset(&(cfg->srcdir), 0, sizeof(float4));
}

void forward_init(Forward* fem, tetmesh* m) {
    fem->mesh = m;
    fem->Ar = NULL;
    fem->Ai = NULL;
    fem->Dr = NULL;
    fem->Di = NULL;
    fem->deldotdel = NULL;
}

void jacobian_init(Jacobian* jac) {
    jac->isnodal = 1;
    jac->nsd = 0;
    jac->nsdcol = 0;
    jac->nn = 0;
    jac->Jmuar = NULL;
    jac->Jmuai = NULL;
    jac->Jdr = NULL;
    jac->Jdi = NULL;
    jac->Phir = NULL;
    jac->Phii = NULL;
    jac->sd = NULL;
    jac->deldotdel = NULL;
    jac->elemid = NULL;
    jac->elembary = NULL;
    jac->relem = NULL;
}

void mesh_init(tetmesh* mesh) {
    mesh->nn = 0;
    mesh->ne = 0;
    mesh->nf = 0;
    mesh->ntype = 0;
    mesh->prop = 0;
    mesh->node = NULL;
    mesh->elem = NULL;
    mesh->face = NULL;
    mesh->type = NULL;
    mesh->med = NULL;
    mesh->evol = NULL;
    mesh->area = NULL;
    mesh->facecenter = NULL;
    mesh->facenormal = NULL;
    mesh->facer = NULL;
    mesh->rows = NULL;
    mesh->cols = NULL;
    mesh->idxcount = NULL;
    mesh->idxsum = NULL;
}

void mesh_clear(tetmesh* mesh) {
    mesh->nn = 0;
    mesh->ne = 0;
    mesh->nf = 0;
    mesh->ntype = 0;
    mesh->ntot = 0;
    mesh->prop = 0;

    if (mesh->node) {
        free(mesh->node);
        mesh->node = NULL;
    }

    if (mesh->elem) {
        free(mesh->elem);
        mesh->elem = NULL;
    }

    if (mesh->face) {
        free(mesh->face);
        mesh->face = NULL;
    }

    if (mesh->area) {
        free(mesh->area);
        mesh->area = NULL;
    }

    if (mesh->type) {
        free(mesh->type);
        mesh->type = NULL;
    }

    if (mesh->med) {
        free(mesh->med);
        mesh->med = NULL;
    }

    if (mesh->evol) {
        free(mesh->evol);
        mesh->evol = NULL;
    }

    if (mesh->facecenter) {
        free(mesh->facecenter);
        mesh->facecenter = NULL;
    }

    if (mesh->facenormal) {
        free(mesh->facenormal);
        mesh->facenormal = NULL;
    }

    if (mesh->facer) {
        free(mesh->facer);
        mesh->facer = NULL;
    }

    if (mesh->rows) {
        free(mesh->rows);
        mesh->rows = NULL;
    }

    if (mesh->cols) {
        free(mesh->cols);
        mesh->cols = NULL;
    }

    if (mesh->idxcount) {
        free(mesh->idxcount);
        mesh->idxcount = NULL;
    }

    if (mesh->idxsum) {
        free(mesh->idxsum);
        mesh->idxsum = NULL;
    }
}

/*
 * Adjoint Jacobian for DOT (mua,D) and MWT (eps_r,sigma).
 *
 * The element/nodal mass-weighted integration here computes
 *    Jmua = -<phi_s, phi_r>_M           (complex when phi is complex)
 *    Jd   = -<grad phi_s, grad phi_r>_K
 * regardless of whether phi solves the diffusion or the Helmholtz
 * equation - the discrete bilinear forms are identical.  In MWT mode
 * the same output represents Jk^2 = -Jmua and the matlab caller chains
 *    Jeps   = -(omega^2*mu0*eps0) * Jmua
 *    Jsigma = ( 1j*omega*mu0    ) * Jmua
 * (see rbmultispectral.m, ishelmholtz branch).  Jd is meaningful only
 * for DOT and should be ignored in MWT mode.
 */
void rb_femjacobian(Config* cfg, tetmesh* mesh, Jacobian* jac) {
    int t, sd, sid, rid, i, j;
    int pairs[3][10] = {{0, 0, 0, 1, 1, 2}, {1, 2, 3, 2, 3, 3}, {1, 2, 3, 5, 6, 8}}; //pairs[0]<->[1], local node pairs, [3] pos in deldotdel
    int inode[4] = {0, 4, 7, 9}; // position in deldotdel for diagonal

    int colnum = (jac->nn ? jac->nn : (jac->isnodal ? mesh->nn : mesh->ne));

    if (jac->Jmuar) {
        memset(jac->Jmuar, 0, jac->nsd * colnum * sizeof(double));
    }

    if (jac->Jmuai) {
        memset(jac->Jmuai, 0, jac->nsd * colnum * sizeof(double));
    }

    if (jac->Jdr) {
        memset(jac->Jdr, 0, jac->nsd * colnum * sizeof(double));
    }

    if (jac->Jdi) {
        memset(jac->Jdi, 0, jac->nsd * colnum * sizeof(double));
    }

    // loop over all elements
    for (t = 0; t < mesh->ne; t++) {
        int* ee = (int*)(mesh->elem + t);

        // first calculate real-part of Jmua and Jd
        for (sd = 0; sd < jac->nsd; sd++) { // loop over sd pairs
            double jdr = 0.0, jmuar = 0.0, dtemp;
            sid = (int)jac->sd[sd] - 1;
            rid = (int)jac->sd[sd + jac->nsd] - 1;

            // accummulate diagonal terms
            for (i = 0; i < 4; i++) {
                dtemp = jac->Phir[sid * mesh->nn + ee[i]] * jac->Phir[rid * mesh->nn + ee[i]];

                if (jac->Phii) {
                    dtemp -= jac->Phii[sid * mesh->nn + ee[i]] * jac->Phii[rid * mesh->nn + ee[i]];
                }

                jmuar += dtemp;
                jdr += jac->deldotdel[inode[i] * mesh->ne + t] * dtemp;
            }

            // accummulate off-diagonal terms
            for (i = 0; i < 6; i++) {
                dtemp = jac->Phir[sid * mesh->nn + ee[pairs[0][i]]] * jac->Phir[rid * mesh->nn + ee[pairs[1][i]]] +
                        jac->Phir[rid * mesh->nn + ee[pairs[0][i]]] * jac->Phir[sid * mesh->nn + ee[pairs[1][i]]];

                if (jac->Phii)
                    dtemp -= jac->Phii[sid * mesh->nn + ee[pairs[0][i]]] * jac->Phii[rid * mesh->nn + ee[pairs[1][i]]] +
                             jac->Phii[rid * mesh->nn + ee[pairs[0][i]]] * jac->Phii[sid * mesh->nn + ee[pairs[1][i]]];

                jmuar += dtemp * 0.5;
                jdr += jac->deldotdel[pairs[2][i] * mesh->ne + t] * dtemp;
            }

            jmuar *= -0.1;
            jdr = -jdr;

            if (jac->nn && jac->nn < mesh->nn) {
                jmuar *= 0.25 * mesh->evol[t];

                for (i = 0; i < 4; i++) {
                    int eid = (int)(jac->elemid[ee[i]] - 0.5);

                    if (jac->elemid[ee[i]] != jac->elemid[ee[i]]) { /*if elemid is nan, the forward node is outside of recon mesh*/
                        continue;
                    }

                    for (j = 0; j < 4; j++) {
                        jac->Jmuar[(int)(jac->relem[j * jac->ne + eid] - 0.5)*jac->nsd + sd] += jmuar * jac->elembary[j * mesh->nn + ee[i]];
                    }
                }
            } else if (jac->isnodal) {
                for (i = 0; i < 4; i++) {
                    jac->Jmuar[ee[i]*jac->nsd + sd] += jmuar * 0.25 * mesh->evol[t];
                }
            } else { // element-based jacobian
                jac->Jmuar[t * jac->nsd + sd] += jmuar * mesh->evol[t];
            }

            if (jac->Jdr) {
                if (jac->nn && jac->nn < mesh->nn) { // support dual-mesh
                    jdr *= 0.25;

                    for (i = 0; i < 4; i++) {
                        int eid = (int)(jac->elemid[ee[i]] - 0.5);

                        if (jac->elemid[ee[i]] != jac->elemid[ee[i]]) { /*if elemid is nan, the forward node is outside of recon mesh*/
                            continue;
                        }

                        for (j = 0; j < 4; j++) {
                            jac->Jdr[(int)(jac->relem[j * jac->ne + eid] - 0.5)*jac->nsd + sd] += jdr * jac->elembary[j * mesh->nn + ee[i]];
                        }
                    }
                } else if (jac->isnodal) {
                    for (i = 0; i < 4; i++) {
                        jac->Jdr[ee[i]*jac->nsd + sd] += jdr * 0.25;
                    }
                } else { // element-based jacobian
                    jac->Jdr[t * jac->nsd + sd] += jdr;
                }
            }
        }
    }

    // if measurement has complex val, then calculate the imag part of J
    if (jac->Phii) {
        for (t = 0; t < mesh->ne; t++) {
            int* ee = (int*)(mesh->elem + t);

            // first calculate real-part of Jmua and Jd
            for (sd = 0; sd < jac->nsd; sd++) { // loop over sd pairs
                double jdi = 0.0, jmuai = 0.0, dtemp;
                sid = (int)jac->sd[sd] - 1;
                rid = (int)jac->sd[sd + jac->nsd] - 1;

                // accummulate diagonal terms
                for (i = 0; i < 4; i++) {
                    dtemp = jac->Phir[sid * mesh->nn + ee[i]] * jac->Phii[rid * mesh->nn + ee[i]];

                    if (jac->Phii) {
                        dtemp += jac->Phii[sid * mesh->nn + ee[i]] * jac->Phir[rid * mesh->nn + ee[i]];
                    }

                    jmuai += dtemp;
                    jdi += jac->deldotdel[inode[i] * mesh->ne + t] * dtemp;
                }

                // accummulate off-diagonal terms
                for (i = 0; i < 6; i++) {
                    dtemp = jac->Phir[sid * mesh->nn + ee[pairs[0][i]]] * jac->Phii[rid * mesh->nn + ee[pairs[1][i]]] +
                            jac->Phir[rid * mesh->nn + ee[pairs[0][i]]] * jac->Phii[sid * mesh->nn + ee[pairs[1][i]]];

                    if (jac->Phii)
                        dtemp += jac->Phii[sid * mesh->nn + ee[pairs[0][i]]] * jac->Phir[rid * mesh->nn + ee[pairs[1][i]]] +
                                 jac->Phii[rid * mesh->nn + ee[pairs[0][i]]] * jac->Phir[sid * mesh->nn + ee[pairs[1][i]]];

                    jmuai += dtemp * 0.5;
                    jdi += jac->deldotdel[pairs[2][i] * mesh->ne + t] * dtemp;
                }

                jmuai *= -0.1;
                jdi = -jdi;

                if (jac->nn && jac->nn < mesh->nn) {
                    jmuai *= 0.25 * mesh->evol[t];

                    for (i = 0; i < 4; i++) {
                        int eid = (int)(jac->elemid[ee[i]] - 0.5);

                        if (jac->elemid[ee[i]] != jac->elemid[ee[i]]) { /*if elemid is nan, the forward node is outside of recon mesh*/
                            continue;
                        }

                        for (j = 0; j < 4; j++) {
                            jac->Jmuai[(int)(jac->relem[j * jac->ne + eid] - 0.5)*jac->nsd + sd] += jmuai * jac->elembary[j * mesh->nn + ee[i]];
                        }
                    }
                } else if (jac->isnodal) { // nodal-based Jacobian
                    for (i = 0; i < 4; i++) {
                        jac->Jmuai[ee[i]*jac->nsd + sd] += jmuai * 0.25 * mesh->evol[t];
                    }
                } else { // element-based jacobian
                    jac->Jmuai[t * jac->nsd + sd] += jmuai * mesh->evol[t];
                }

                if (jac->Jdi) {
                    if (jac->nn && jac->nn < mesh->nn) { // support dual-mesh
                        jdi *= 0.25;

                        for (i = 0; i < 4; i++) {
                            int eid = (int)(jac->elemid[ee[i]] - 0.5);

                            if (jac->elemid[ee[i]] != jac->elemid[ee[i]]) { /*if elemid is nan, the forward node is outside of recon mesh*/
                                continue;
                            }

                            for (j = 0; j < 4; j++) {
                                jac->Jdi[(int)(jac->relem[j * jac->ne + eid] - 0.5)*jac->nsd + sd] += jdi * jac->elembary[j * mesh->nn + ee[i]];
                            }
                        }
                    } else if (jac->isnodal) {
                        for (i = 0; i < 4; i++) {
                            jac->Jdi[ee[i]*jac->nsd + sd] += jdi * 0.25;
                        }
                    } else { // element-based jacobian
                        jac->Jdi[t * jac->nsd + sd] += jdi;
                    }
                }
            }
        }
    }
}

void rb_deldotdel(Config* cfg, tetmesh* mesh, Forward* fem) {
    float3* N;
    int i1, i2, i3, i4, i, j, t;
    double derx[4], dery[4], derz[4];
    double Ve, RVe6;

    /* map coordinate pointers*/
    N = mesh->node;

    /* loop over all elements*/
    for (t = 0; t < mesh->ne; t++) {
        RVe6 = 1.0 / (mesh->evol[t] * 6.0);
        Ve = mesh->evol[t];

        i1 = mesh->elem[t].x;
        i2 = mesh->elem[t].y;
        i3 = mesh->elem[t].z;
        i4 = mesh->elem[t].w;

        /*calculate del_x,del_y and del_z for 4 basis functions per elem*/
        derx[0] = -((N[i3].y * N[i4].z - N[i3].z * N[i4].y) - N[i2].y * (N[i4].z - N[i3].z)
                    + N[i2].z * (N[i4].y - N[i3].y)) * RVe6;
        dery[0] = ((N[i3].x * N[i4].z - N[i4].x * N[i3].z) - N[i2].x * (N[i4].z - N[i3].z)
                   + N[i2].z * (N[i4].x - N[i3].x)) * RVe6;
        derz[0] = -((N[i3].x * N[i4].y - N[i3].y * N[i4].x) - N[i2].x * (N[i4].y - N[i3].y)
                    + N[i2].y * (N[i4].x - N[i3].x)) * RVe6;

        derx[1] = ((N[i3].y * N[i4].z - N[i3].z * N[i4].y) - N[i1].y * (N[i4].z - N[i3].z)
                   + N[i1].z * (N[i4].y - N[i3].y)) * RVe6;
        dery[1] = -((N[i3].x * N[i4].z - N[i4].x * N[i3].z) - N[i1].x * (N[i4].z - N[i3].z)
                    + N[i1].z * (N[i4].x - N[i3].x)) * RVe6;
        derz[1] = ((N[i3].x * N[i4].y - N[i3].y * N[i4].x) - N[i1].x * (N[i4].y - N[i3].y)
                   + N[i1].y * (N[i4].x - N[i3].x)) * RVe6;

        derx[2] = -((N[i2].y * N[i4].z - N[i2].z * N[i4].y) - N[i1].y * (N[i4].z - N[i2].z)
                    + N[i1].z * (N[i4].y - N[i2].y)) * RVe6;
        dery[2] = ((N[i2].x * N[i4].z - N[i4].x * N[i2].z) - N[i1].x * (N[i4].z - N[i2].z)
                   + N[i1].z * (N[i4].x - N[i2].x)) * RVe6;
        derz[2] = -((N[i2].x * N[i4].y - N[i2].y * N[i4].x) - N[i1].x * (N[i4].y - N[i2].y)
                    + N[i1].y * (N[i4].x - N[i2].x)) * RVe6;

        derx[3] = ((N[i2].y * N[i3].z - N[i2].z * N[i3].y) - N[i1].y * (N[i3].z - N[i2].z)
                   + N[i1].z * (N[i3].y - N[i2].y)) * RVe6;
        dery[3] = -((N[i2].x * N[i3].z - N[i3].x * N[i2].z) - N[i1].x * (N[i3].z - N[i2].z)
                    + N[i1].z * (N[i3].x - N[i2].x)) * RVe6;
        derz[3] = ((N[i2].x * N[i3].y - N[i2].y * N[i3].x) - N[i1].x * (N[i3].y - N[i2].y)
                   + N[i1].y * (N[i3].x - N[i2].x)) * RVe6;

        /*loop over upper triangle*/
        for (i = 0; i < 4; i++)
            for (j = i; j < 4; j++) {
                fem->deldotdel[t * 10 + (i << 2) + j - i - (i > 1) - ((i > 2) << 1)] = (derx[i] * derx[j] + dery[i] * dery[j] + derz[i] * derz[j]) * Ve;
            }
    }
}

void rb_femmatrix_elem(Config* cfg, tetmesh* mesh, Forward* fem) {
    double mua, musp, dcoeff, nref;
    int i, j, t, ii, jj, ij, *ee;
    double ra, ri, Ve, deldotdel, sm;

    /*initialize output variables*/
    if (fem->Dr) {
        memset(fem->Dr, 0, mesh->nn * sizeof(double));
    }

    if (fem->Di) {
        memset(fem->Di, 0, mesh->nn * sizeof(double));
    }

    if (fem->Ar) {
        memset(fem->Ar, 0, mesh->ntot * sizeof(double));
    }

    if (fem->Ai) {
        memset(fem->Ai, 0, mesh->ntot * sizeof(double));
    }

    /* loop over all elements*/
    for (t = 0; t < mesh->ne; t++) {
        Ve = mesh->evol[t];
        ee = (int*)(mesh->elem + t);

        mua = mesh->med[mesh->type[t]].mua;
        musp = mesh->med[mesh->type[t]].mus;
        nref = mesh->med[mesh->type[t]].n;

        dcoeff = 1.0 / (3.0 * (mua + musp));

        /*loop over index i*/
        for (i = 0; i < 4; i++) {
            ii = ee[i];

            /*loop over index j, upper triangle only*/
            for (j = i; j < 4; j++) {
                jj = ee[j];
                deldotdel = fem->deldotdel[t * 10 + (i << 2) + j - i - (i > 1) - ((i > 2) << 1)];

                sm = Ve * 0.05;

                if (i == j) {
                    sm = Ve * 0.1;
                }

                /* D*Ve*delFi*delF */
                ra = dcoeff * deldotdel + mua * sm;
                ri = cfg->omega * R_C0 * nref * sm;

                /*add the values to the vectorized sparse matrix*/
                if (ii == jj) {
                    fem->Dr[ii] += ra;

                    if (fem->Di) {
                        fem->Di[ii] += ri;
                    }
                } else {
                    ij = sub2ind(mesh, ii, jj);
                    fem->Ar[ij] += ra;

                    if (fem->Ai) {
                        fem->Ai[ij] += ri;
                    }
                }
            }
        }
    }
}

void rb_femmatrix_nodal(Config* cfg, tetmesh* mesh, Forward* fem) {
    double mua, musp, dcoeff, nref;
    int i, j, k, t, ii, jj, ij, *ee;
    double ra, ri, Ve, deldotdel, sm;

    /*initialize output variables*/
    if (fem->Dr) {
        memset(fem->Dr, 0, mesh->nn * sizeof(double));
    }

    if (fem->Di) {
        memset(fem->Di, 0, mesh->nn * sizeof(double));
    }

    if (fem->Ar) {
        memset(fem->Ar, 0, mesh->ntot * sizeof(double));
    }

    if (fem->Ai) {
        memset(fem->Ai, 0, mesh->ntot * sizeof(double));
    }

    /* loop over all elements*/
    for (t = 0; t < mesh->ne; t++) {
        Ve = mesh->evol[t];
        ee = (int*)(mesh->elem + t);

        /*loop over index i*/
        for (i = 0; i < 4; i++) {
            ii = ee[i];

            /*loop over index j*/
            for (j = i; j < 4; j++) {
                jj = ee[j];
                deldotdel = fem->deldotdel[t * 10 + (i << 2) + j - i - (i > 1) - ((i > 2) << 1)] * 0.25;

                ra = 0.0;
                ri = 0.0;

                for (k = 0; k < 4; k++) {
                    sm = Ve * (1.0 / 120.);

                    if (i == j || i == k || j == k) {
                        sm = Ve * (1.0 / 60.);
                    }

                    if (i == j && i == k) {
                        sm = Ve * (1.0 / 20.);
                    }

                    mua = mesh->med[mesh->type[ee[k]]].mua;
                    musp = mesh->med[mesh->type[ee[k]]].mus;
                    nref = mesh->med[mesh->type[ee[k]]].n;

                    dcoeff = 1.0 / (3.0 * (mua + musp));

                    /* D*Ve*delFi*delF */
                    ra += dcoeff * deldotdel + mua * sm;
                    ri += cfg->omega * R_C0 * nref * sm;
                }

                /*add the values to the vectorized sparse matrix*/
                if (ii == jj) {
                    fem->Dr[ii] += ra;

                    if (fem->Di) {
                        fem->Di[ii] += ri;
                    }
                } else {
                    ij = sub2ind(mesh, ii, jj);
                    fem->Ar[ij] += ra;

                    if (fem->Ai) {
                        fem->Ai[ij] += ri;
                    }
                }
            }
        }
    }
}


void rb_fem_bc(Config* cfg, tetmesh* mesh, Forward* fem) {
    int t;
    double val;
    double Reff = (1.0 - cfg->reff) / (1.0 + cfg->reff) * (1.0 / 12.0);

    for (t = 0; t < mesh->nf; t++) {
        val = mesh->area[t] * Reff;
        fem->Dr[mesh->face[t].x] += val;
        fem->Dr[mesh->face[t].y] += val;
        fem->Dr[mesh->face[t].z] += val;

        val = val * 0.5;
        fem->Ar[sub2ind(mesh, mesh->face[t].x, mesh->face[t].y)] += val;
        fem->Ar[sub2ind(mesh, mesh->face[t].x, mesh->face[t].z)] += val;
        fem->Ar[sub2ind(mesh, mesh->face[t].y, mesh->face[t].z)] += val;
    }
}

/*
 * Element-based Helmholtz/MWT stiffness matrix.
 * Per-element coefficients (from rbfemlhs.m, ishelmholtz branch):
 *   avol  = 1                                        (stiffness)
 *   breal = -omega^2 * mu0 * eps0 * eps_r            (mass, real part)
 *   bimag =  omega   * mu0 * sigma                   (mass, imag part)
 * Property slots are reused: mua->eps_r, mus->sigma, g->mu0, n->n.
 * The mass-matrix weights (0.10*Ve diag, 0.05*Ve off-diag) and sparsity
 * are identical to the DOT element path, so the routine mirrors
 * rb_femmatrix_elem with the (avol,breal,bimag) substitution.
 */
void rb_femmatrix_helmholtz_elem(Config* cfg, tetmesh* mesh, Forward* fem) {
    double eps_r, sigma, mu0_e;
    double breal, bimag;
    int i, j, t, ii, jj, ij, *ee;
    double ra, ri, Ve, deldotdel, sm;
    double omega = cfg->omega;
    double omega2 = omega * omega;

    if (fem->Dr) {
        memset(fem->Dr, 0, mesh->nn * sizeof(double));
    }

    if (fem->Di) {
        memset(fem->Di, 0, mesh->nn * sizeof(double));
    }

    if (fem->Ar) {
        memset(fem->Ar, 0, mesh->ntot * sizeof(double));
    }

    if (fem->Ai) {
        memset(fem->Ai, 0, mesh->ntot * sizeof(double));
    }

    for (t = 0; t < mesh->ne; t++) {
        Ve = mesh->evol[t];
        ee = (int*)(mesh->elem + t);

        eps_r = mesh->med[mesh->type[t]].mua;
        sigma = mesh->med[mesh->type[t]].mus;
        mu0_e = mesh->med[mesh->type[t]].g;

        breal = -omega2 * mu0_e * EPS0_MM * eps_r;
        bimag =  omega  * mu0_e * sigma;

        for (i = 0; i < 4; i++) {
            ii = ee[i];

            for (j = i; j < 4; j++) {
                jj = ee[j];
                deldotdel = fem->deldotdel[t * 10 + (i << 2) + j - i - (i > 1) - ((i > 2) << 1)];

                sm = Ve * 0.05;

                if (i == j) {
                    sm = Ve * 0.1;
                }

                /* stiffness (avol=1) + real part of -k^2 mass */
                ra = deldotdel + breal * sm;
                /* imag part of -k^2 mass */
                ri = bimag * sm;

                if (ii == jj) {
                    fem->Dr[ii] += ra;

                    if (fem->Di) {
                        fem->Di[ii] += ri;
                    }
                } else {
                    ij = sub2ind(mesh, ii, jj);
                    fem->Ar[ij] += ra;

                    if (fem->Ai) {
                        fem->Ai[ij] += ri;
                    }
                }
            }
        }
    }
}

/*
 * Nodal Helmholtz/MWT stiffness matrix.  Same sparsity & mass weights as
 * the DOT nodal path; only the per-node (avol,breal,bimag) substitution
 * differs.  Stiffness contribution per element is summed over its 4
 * nodes with weight 0.25 (already folded into deldotdel here), so with
 * avol=1 the four-fold sum recovers the full deldotdel.
 */
void rb_femmatrix_helmholtz_nodal(Config* cfg, tetmesh* mesh, Forward* fem) {
    double eps_r, sigma, mu0_e;
    double breal_k, bimag_k;
    int i, j, k, t, ii, jj, ij, *ee;
    double ra, ri, Ve, deldotdel, sm;
    double omega = cfg->omega;
    double omega2 = omega * omega;

    if (fem->Dr) {
        memset(fem->Dr, 0, mesh->nn * sizeof(double));
    }

    if (fem->Di) {
        memset(fem->Di, 0, mesh->nn * sizeof(double));
    }

    if (fem->Ar) {
        memset(fem->Ar, 0, mesh->ntot * sizeof(double));
    }

    if (fem->Ai) {
        memset(fem->Ai, 0, mesh->ntot * sizeof(double));
    }

    for (t = 0; t < mesh->ne; t++) {
        Ve = mesh->evol[t];
        ee = (int*)(mesh->elem + t);

        for (i = 0; i < 4; i++) {
            ii = ee[i];

            for (j = i; j < 4; j++) {
                jj = ee[j];
                deldotdel = fem->deldotdel[t * 10 + (i << 2) + j - i - (i > 1) - ((i > 2) << 1)] * 0.25;

                ra = 0.0;
                ri = 0.0;

                for (k = 0; k < 4; k++) {
                    sm = Ve * (1.0 / 120.);

                    if (i == j || i == k || j == k) {
                        sm = Ve * (1.0 / 60.);
                    }

                    if (i == j && i == k) {
                        sm = Ve * (1.0 / 20.);
                    }

                    eps_r = mesh->med[mesh->type[ee[k]]].mua;
                    sigma = mesh->med[mesh->type[ee[k]]].mus;
                    mu0_e = mesh->med[mesh->type[ee[k]]].g;

                    breal_k = -omega2 * mu0_e * EPS0_MM * eps_r;
                    bimag_k =  omega  * mu0_e * sigma;

                    /* avol=1 stiffness (4*0.25*deldotdel = deldotdel) + nodal mass */
                    ra += deldotdel + breal_k * sm;
                    ri += bimag_k * sm;
                }

                if (ii == jj) {
                    fem->Dr[ii] += ra;

                    if (fem->Di) {
                        fem->Di[ii] += ri;
                    }
                } else {
                    ij = sub2ind(mesh, ii, jj);
                    fem->Ar[ij] += ra;

                    if (fem->Ai) {
                        fem->Ai[ij] += ri;
                    }
                }
            }
        }
    }
}

/*
 * First-order Bayliss-Turkel radiation boundary condition for MWT/
 * Helmholtz (matches rbfemlhs.m).  For each boundary face t the per-
 * face contribution is
 *   bccoef(t) = (i*kbg - 1/(2*r(t))) * (rhat dot n_hat)
 * where kbg = sqrt(omega^2*mu*eps0*eps_bulk - i*omega*mu*sigma_bulk)
 * is the bulk wavenumber and r(t)=|facecenter(t)-rbcorigin|.
 * The face mass-matrix integrates to (area/6) on the three diagonal
 * entries and (area/12)=(area/6)*0.5 on the three off-diagonal pairs.
 */
void rb_fem_bc_helmholtz(Config* cfg, tetmesh* mesh, Forward* fem) {
    if (!mesh->facecenter || !mesh->facenormal || !mesh->facer) {
        MEXERROR("Helmholtz BC requires cfg.facecenter, cfg.facenormal and cfg.facer");
    }

    typedef std::complex<double> cplx;

    double omega = cfg->omega;
    double eps_b = cfg->bulkprop.mua;
    double sigma_b = cfg->bulkprop.mus;
    double mu_b = cfg->bulkprop.g;

    /* k^2 of the bulk medium and kbg = sqrt(k^2). std::sqrt picks the
     * principal branch (Re kbg >= 0), matching matlab's sqrt(). */
    cplx k2bg(omega * omega * mu_b * EPS0_MM * eps_b,
              -omega * mu_b * sigma_b);
    cplx kbg = std::sqrt(k2bg);

    for (int t = 0; t < mesh->nf; t++) {
        double rx = mesh->facecenter[t].x - cfg->rbcorigin.x;
        double ry = mesh->facecenter[t].y - cfg->rbcorigin.y;
        double rz = mesh->facecenter[t].z - cfg->rbcorigin.z;
        double r = mesh->facer[t];

        if (r == 0.0) {
            continue;    /* face centroid coincides with origin: skip */
        }

        double rdotn = (rx * mesh->facenormal[t].x
                        + ry * mesh->facenormal[t].y
                        + rz * mesh->facenormal[t].z) / r;

        cplx bccoef = (cplx(0.0, 1.0) * kbg - cplx(1.0 / (2.0 * r), 0.0)) * rdotn;
        cplx adiag = (mesh->area[t] / 6.0) * bccoef;
        cplx aoffd = adiag * 0.5;

        double dr = adiag.real(), di = adiag.imag();
        double or_ = aoffd.real(), oi = aoffd.imag();

        int v0 = mesh->face[t].x;
        int v1 = mesh->face[t].y;
        int v2 = mesh->face[t].z;

        fem->Dr[v0] += dr;
        fem->Dr[v1] += dr;
        fem->Dr[v2] += dr;

        if (fem->Di) {
            fem->Di[v0] += di;
            fem->Di[v1] += di;
            fem->Di[v2] += di;
        }

        int ij;
        ij = sub2ind(mesh, v0, v1);
        fem->Ar[ij] += or_;

        if (fem->Ai) {
            fem->Ai[ij] += oi;
        }

        ij = sub2ind(mesh, v0, v2);
        fem->Ar[ij] += or_;

        if (fem->Ai) {
            fem->Ai[ij] += oi;
        }

        ij = sub2ind(mesh, v1, v2);
        fem->Ar[ij] += or_;

        if (fem->Ai) {
            fem->Ai[ij] += oi;
        }
    }
}

int sub2ind(tetmesh* mesh, int i, int j) {
    int base, ij = -1, k;
    base = (int)(mesh->idxsum[i] - mesh->idxcount[i]);

    for (k = 0; k < (int)(mesh->idxcount[i]); k++) {
        ij = k + base;

        if ((int)(mesh->cols[ij]) == j) {
            break;
        }
    }

    return ij;
}

extern "C" int rb3_throw_exception(const int id, const char* msg, const char* filename, const int linenum) {
    printf("Redbird ERROR (%d): %s in unit %s:%d\n", id, msg, filename, linenum);
    throw (msg);
    return id;
}

void rb3_usage() {
    printf("Usage:\n"
           "    [Adiag,Aoff,deldotdel] = rbfemmatrix(cfg);\n"
           "    [Jmua,Jd]              = rbfemmatrix(cfg,sd,phi[,deldotdel[,isnodal]]);\n"
           "\nDOT mode (default): cfg.prop columns are [mua mus g n].\n"
           "MWT/Helmholtz mode (engaged by cfg.helmholtz=1 or by providing\n"
           "cfg.bulkprop): cfg.prop columns are reinterpreted as\n"
           "[eps_r sigma(S/mm) mu0(H/mm) n].  In MWT mode the caller must\n"
           "also provide cfg.bulkprop=[eps_r sigma mu0 (n)], cfg.rbcorigin,\n"
           "cfg.facecenter, cfg.facenormal and cfg.facer for the first-order\n"
           "Bayliss-Turkel radiation boundary condition.  The Jacobian output\n"
           "Jmua in MWT mode represents -<phi_s,phi_r>_M; the caller chains it\n"
           "to Jeps = -(omega^2*mu0*eps0)*Jmua and Jsigma = (1j*omega*mu0)*Jmua.\n");
}

extern "C" void rb3_flush() {
#ifndef MATLAB_MEX_FILE
    mexEvalString("fflush(stdout);");
#else
    mexEvalString("pause(.0001);");
#endif
}
