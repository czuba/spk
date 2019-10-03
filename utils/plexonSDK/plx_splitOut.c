#include <math.h>
#include "mex.h"
#include <stdio.h>
/* 12-aug-2009*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *filename,*filename1,*filename2;
    FILE *fid, *fid1, *fid2;
    int buflen, buflen1, buflen2;
    short type, upperbyte,channel,unit,nwf, nwords,wf[64];
    int header[64],ndsp,nevents,nslow,timestamp;
    unsigned long en;
    double *corr;
    unsigned long c;
    char inter;
    int i;
 /* Get filename */
    
    buflen = mxGetM(prhs[0]) * mxGetN(prhs[0])+1;
    filename = mxCalloc(buflen, sizeof(char));
    mxGetString(prhs[0], filename, buflen);
 
    /* Outputs */
    if (nrhs<4) {
        /* Make generic output names */
        filename1 = mxCalloc(buflen+2, sizeof(char));
        filename2 = mxCalloc(buflen+2, sizeof(char));
        for(i = 0; i<buflen; i++)
        {
            filename1[i+2          ] = filename[i];
            filename2[i+2          ] = filename[i];
        }
        filename1[1] = 45;
        filename2[1] = 45;
        
        filename1[0] = 97;
        filename2[0] = 98;

    } else {
        /* Use inputs*/
        buflen1 = mxGetM(prhs[2]) * mxGetN(prhs[2])+1;
        filename1 = mxCalloc(buflen1, sizeof(char));
        mxGetString(prhs[2], filename1, buflen1);

        buflen2 = mxGetM(prhs[3]) * mxGetN(prhs[3])+1;
        filename2 = mxCalloc(buflen2, sizeof(char));
        mxGetString(prhs[3], filename2, buflen2);

    }
    
    fid = fopen(filename,"rb");
    if (fid == NULL)
    {
        mexPrintf("PLX File not found\n");
        return;
    }
    
    
    
    fread(&header,4,64,fid);
    ndsp = header[35];
    nevents = header[36];
    nslow = header[37];
    
    fseek(fid,0,SEEK_END);
    en = ftell(fid);
   
    fseek(fid,0,SEEK_SET);
    
    
    fid1 = fopen(filename1,"wb+");
    fid2 = fopen(filename2,"wb+");
    if (fid1 == NULL || fid2 == NULL)
    {
        mexPrintf("Error: Unable to write a new file!\n");
        return;
    }
    
    for(i=0; i<1020*ndsp + 296*nevents + 296*nslow+7504;i++)
    {
        fread (&inter, 1, 1, fid );
        fwrite(&inter, 1, 1, fid1);
        fwrite(&inter, 1, 1, fid2);
    }
    
    corr = mxGetPr(prhs[1]);
    c = corr[0]*10000;
    while ((ftell(fid)!=en || en == -1) && !feof(fid))
    {
        fread(&type,2,1,fid);
        if (feof(fid))
            break;
        fread(&upperbyte,2,1,fid);
        fread(&timestamp,4,1,fid);
        fread(&channel,2,1,fid);
        fread(&unit,2,1,fid);
        fread(&nwf,2,1,fid);
        fread(&nwords,2,1,fid);
        if (nwords > 0)
            fread(&wf,2,nwords,fid);
        if (timestamp<c)
        {
            fwrite(&type,2,1,fid1);
            fwrite(&upperbyte,2,1,fid1);
            fwrite(&timestamp,4,1,fid1);
            fwrite(&channel,2,1,fid1);
            fwrite(&unit,2,1,fid1);
            fwrite(&nwf,2,1,fid1);
            fwrite(&nwords,2,1,fid1);
            if (nwords > 0)
                fwrite(&wf,2,nwords,fid1);
                
        }
        else
        {
            timestamp = timestamp-c;
            fwrite(&type,2,1,fid2);
            fwrite(&upperbyte,2,1,fid2);
            fwrite(&timestamp,4,1,fid2);
            fwrite(&channel,2,1,fid2);
            fwrite(&unit,2,1,fid2);
            fwrite(&nwf,2,1,fid2);
            fwrite(&nwords,2,1,fid2);
            if (nwords > 0)
                fwrite(&wf,2,nwords,fid2);
        } 
    }
    fclose(fid );
    fclose(fid1);
    fclose(fid2);
    
}
