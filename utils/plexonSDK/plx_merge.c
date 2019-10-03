#include <math.h>
#include "mex.h"
// 2012-10-15 TBC: Got it from Amin & added an optional output filename & some matlab outputs (combined filename & split point)
// 2013-03-13 TBC: fixed up with Amin's improvements to skip plx files dumping LFP data in the middle of the spike data.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *filename,*filename1,*filename2;
    FILE *fid, *fid1, *fid2;
    int buflen1,buflen2, buflen0;
    unsigned int corr;
    int i;
    short type, upperbyte,channel,unit,nwf, nwords,wf[500];
    int header[64],ndsp,nevents,nslow,timestamp,ts;
    double *outts;
    char inter;
    
    // Get filename
    
    buflen1 = mxGetM(prhs[0]) * mxGetN(prhs[0])+1;
    filename1 = mxCalloc(buflen1, sizeof(char));
    mxGetString(prhs[0], filename1, buflen1);
    
    buflen2 = mxGetM(prhs[1]) * mxGetN(prhs[1])+1;
    filename2 = mxCalloc(buflen2, sizeof(char));
    mxGetString(prhs[1], filename2, buflen2);
    
    // Output filename
    if (nrhs< 3) {
        mexPrintf("**************\nNo output name provided...creating one from input filenames.\n**************\n");
        
        filename = mxCalloc(buflen1+buflen2-4, sizeof(char));
        for(i = 0; i<buflen1-5; i++)
            filename[i          ] = filename1[i];
        filename[buflen1-5] = 43;
        
        for(i = 0; i<buflen2; i++)
            filename[i+buflen1-4] = filename2[i];
    } else {
        buflen0 = mxGetM(prhs[2]) * mxGetN(prhs[2])+1;
        filename = mxCalloc(buflen0, sizeof(char));
        mxGetString(prhs[2], filename, buflen0);
    }
    
    
    fid1 = fopen(filename1,"rb");
    if (fid1 == NULL)
    {
        mexPrintf("PLX File not found\n");
        return;
    }
    
    fid = fopen(filename,"wb+");
    if (fid == NULL)
    {
        mexPrintf("Unable to open the destination file\n");
        return;
    }
    
    fread(&header,4,64,fid1);
    ndsp = header[35];
    nevents = header[36];
    nslow = header[37];
    
    fseek(fid1,0,SEEK_SET);
    
    for(i=0; i<1020*ndsp + 296*nevents + 296*nslow+7504;i++)
    {
        fread (&inter, 1, 1, fid1);
        fwrite(&inter, 1, 1, fid );
    }
    
    mexPrintf("Reading ");
    mexPrintf(filename1);
    mexPrintf("\n");
    
    while (!feof(fid1))
    {
        fread (&type,2,1,fid1);
        fread (&upperbyte,2,1,fid1);
        fread (&timestamp,4,1,fid1);
        fread (&channel,2,1,fid1);
        fread (&unit,2,1,fid1);
        fread (&nwf,2,1,fid1);
        if (feof(fid1))
            break;
        fread (&nwords,2,1,fid1);
        if (nwords > 0)
        {
            //fread (&wf,2,nwords,fid1);
            if (type !=5)
            {
                fread (&wf,2,nwords,fid1);
            }
            else
            {
                fseek (fid1,2*nwords,SEEK_CUR);
            }
            
        }
        
        
        if (type !=5)
        {
            fwrite(&type,2,1,fid );
            fwrite(&upperbyte,2,1,fid );
            fwrite(&timestamp,4,1,fid );
            fwrite(&channel,2,1,fid );
            fwrite(&unit,2,1,fid );
            fwrite(&nwf,2,1,fid );
            fwrite(&nwords,2,1,fid );
            if (nwords > 0)
            {
                fwrite(&wf,2,nwords,fid );
            }
        }
        
        if (type == 1)
            ts = timestamp;
    }
    
    fclose(fid1);
    
    fid2 = fopen(filename2,"rb");
    if (fid2 == NULL)
    {
        mexPrintf("PLX File not found\n");
        return;
    }
    
    fread(&header,4,64,fid2);
    ndsp = header[35];
    nevents = header[36];
    nslow = header[37];
    
    fseek(fid2,1020*ndsp + 296*nevents + 296*nslow+7504,SEEK_SET);
    
    corr = 4+(long)ts/10000;
    mexPrintf("Reading ");
    mexPrintf(filename2);
    mexPrintf("\n");
    
    while (!feof(fid2))
    {
        fread (&type,2,1,fid2);
        fread (&upperbyte,2,1,fid2);
        fread (&timestamp,4,1,fid2);
        fread (&channel,2,1,fid2);
        fread (&unit,2,1,fid2);
        fread (&nwf,2,1,fid2);
        if (feof(fid2))
            break;
        fread (&nwords,2,1,fid2);
        if (nwords > 0)
        {
            //fread (&wf,2,nwords,fid2);
            if (type !=5)
            {
                fread (&wf,2,nwords,fid2);
            }
            else
            {
                fseek (fid2,2*nwords,SEEK_CUR);
            }
        }
        
        if (type !=5)
        {
            fwrite(&type,2,1,fid );
            fwrite(&upperbyte,2,1,fid );
            timestamp += corr*10000;
            fwrite(&timestamp,4,1,fid );
            fwrite(&channel,2,1,fid );
            fwrite(&unit,2,1,fid );
            fwrite(&nwf,2,1,fid );
            fwrite(&nwords,2,1,fid );
            if (nwords > 0)
            {
                fwrite(&wf,2,nwords,fid );
            }
        }
    }
    
    fclose(fid2);
    fclose(fid );
    mexPrintf("Saved output as: %s\n",filename);
    mexPrintf("mergeindex: %1d\n",corr);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    outts = mxGetPr(plhs[0]);
    outts[0] = (double)corr;
    
}