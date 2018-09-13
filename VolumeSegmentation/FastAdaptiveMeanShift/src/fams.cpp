/////////////////////////////////////////////////////////////////////////////
// Name:        fams.cpp
// Purpose:     fast adaptive meanshift implementation
// Author:      Ilan Shimshoni
// Modified by: Bogdan Georgescu
// Created:     08/14/2003
// Version:     v0.1
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include "fams.h"


FAMS::FAMS(int no_lsh)
{
   hasPoints_ = 0;
   nsel_ = 0;
   npm_ = 0;
   hashCoeffs_ = NULL;
   noLSH_=no_lsh;

   time_t tt1;
   time(&tt1);
   srand((unsigned)tt1);
//   srand(100);

   nnres1_ = 0;
   nnres2_ = 0;
}

FAMS::~FAMS()
{
   CleanData();
   
}

void FAMS::CleanData()
{
   CleanPoints();
   CleanSelected();
   CleanPrunedModes();
   CleanHash();
}

void FAMS::CleanPoints()
{
   if (hasPoints_)
   {
      delete [] points_;
      delete [] data_;
      delete [] rr_;
      hasPoints_ = 0;
   }
}

void FAMS::CleanSelected()
{
   if (nsel_ > 0)
   {
      delete [] psel_;
      delete [] modes_;
      delete [] hmodes_;
      nsel_ = 0;
   }
}

void FAMS::CleanPrunedModes()
{
   if (npm_ > 0)
   {
      delete [] prunedmodes_;
      delete [] nprunedmodes_;
      npm_ = 0;
   }
}

void FAMS::CleanHash()
{
   if (hashCoeffs_ != NULL)
   {
      delete [] hashCoeffs_;
      hashCoeffs_ = NULL;
   }
}

int FAMS::LoadPoints(char* filename)
{
   bgLog("Load data points from %s...",filename);
   CleanPoints();
   FILE* fd;
   fd = fopen(filename, "r");

   if (!fd)
   {
      bgLog("Error opening %s\n", filename);
      return 1;
   }
   fscanf(fd, "%d %d", &n_, &d_);
   if ((n_<1) || (d_<1))
   {
      bgLog("Error reading %s\n", filename);
      fclose(fd);
      return 1;
   }

   // allocate data
   float* pttemp;
   pttemp = new float[n_*d_];
   int i;
   for (i=0; (i<(n_*d_)) && (fscanf(fd, "%g", &pttemp[i]) == 1); i++);
   fclose(fd);

   if (i!= (n_*d_))
   {
      bgLog("Error reading %s\n", filename);
      delete [] pttemp;
      return 1;
   }

   // allocate and convert to integer
   for (i=0, minVal_=pttemp[0], maxVal_=pttemp[0]; i<(n_*d_); i++)
   {
      if (minVal_>pttemp[i])
         minVal_ = pttemp[i];
      else if (maxVal_<pttemp[i])
         maxVal_ = pttemp[i];
   }
   data_ = new unsigned short[n_*d_];
   rr_ = new double[d_];
   hasPoints_ = 1;
   float deltaVal = maxVal_-minVal_;
   if (deltaVal == 0) deltaVal = 1;
   for (i=0; i<(n_*d_); i++)
      data_[i] = (unsigned short) (65535.0*(pttemp[i]-minVal_)/deltaVal);
   delete [] pttemp;
   dataSize_ = d_*sizeof(unsigned short);

   points_ = new fams_point[n_];
   unsigned short* dtempp;
   for (i=0, dtempp=data_; i<n_; i++, dtempp+=d_)
   {
      points_[i].data_ = dtempp;
      points_[i].usedFlag_ = 0;
   }
   bgLog("done\n");
   return 0;
}


// Choose a subset of points on which to perform the mean shift operation 
void FAMS::SelectMsPoints(double percent,int jump)
{
   if (hasPoints_ == 0)
      return;
   int i;
   int tsel;
   if(percent){
      tsel = (int)(n_*percent/100.0);
      if (tsel != nsel_)
      {
         CleanSelected();
         nsel_ = tsel;
         psel_ = new int[nsel_];
         modes_ = new unsigned short[nsel_*d_];
         hmodes_ = new unsigned int[nsel_];
      }
      for(i=0; i<nsel_;  i++)
         psel_[i] = MyMin(n_-1,(int)(drand48()*n_));
   }
   else
   {
      tsel = (int)ceil(n_/(jump+0.0));
      if (tsel != nsel_)
      {
         CleanSelected();
         nsel_ = tsel;
         psel_ = new int[nsel_];
         modes_ = new unsigned short[nsel_*d_];
         hmodes_ = new unsigned int[nsel_];
      }
      for(i=0; i<nsel_; i++)
         psel_[i]=i*jump;
   }
}

void FAMS::MakeCuts(fams_partition* cuts)
{
   int i;
   for(i=0; i<L_; i++)
      MakeCutL(cuts[i]);
}

// data-driven uniform partition

void FAMS::MakeCutL(fams_partition& cut)
{
   int n1 = (int)floor(K_/(1.0*d_));
   int i,j;
   int ncu=0;
   int w;
   for(i=0; i<d_; i++)
   {
      for(j=0; j<n1; j++)
      {
         cut[ncu].which_ = i;
         w = MyMin((int)(drand48()*n_),n_-1);
         cut[ncu].where_ = points_[w].data_[i];
         ncu++;
      }
   }
   /*
   int *which = new int[d_];
   for (i=0; i<d_; i++)
      which[i]=i;
   int wh, ndleft, itmp;
   ndleft = d_;
   for (i=0; i<(K_-ncu); i++)
   {
      wh = MyMin((int)(drand48()*ndleft),ndleft-1);
      itmp = which[wh+i];
      which[wh+i] = which[i];
      which[i] = itmp;
      ndleft--;
   }

   for(i=0 ; ncu < K_; ncu++, i++)
   {
      w = MyMin((int)(drand48()*n_),n_-1);
      cut[ncu].which_ = which[i];
      cut[ncu].where_ = points_[w].data_[which[i]];
   }
   delete [] which;
   */
   int *which = new int[d_];
   memset(which,0,sizeof(int)*d_);
   for (; ncu<K_; )
   {
      int wh = MyMin((int)(drand48()*d_),d_-1);
      if(which[wh]!=0)
         continue;
      which[wh]=1;
      int w = MyMin((int)(drand48()*n_), n_-1);
      cut[ncu].which_ = wh;
      cut[ncu].where_ = points_[w].data_[wh];
      ncu++;
   }
   delete [] which;
}

void FAMS::InitHash(int nk)
{
   CleanHash();
   hashCoeffs_ = new int[nk];
   for(int i=0; i<nk; i++)
      hashCoeffs_[i] = rand();
}

/* compute the hash key and and the double hash key if needed 
It is possible to give M the the size of the hash table and 
hjump for the double hash key 
*/

int FAMS::HashFunction(int *cutVals, int whichPartition, int kk,int M,int *hjump)
{
   int i;
   int res = whichPartition;
   for(i=0; i<kk; i++)
   {
      res += cutVals[i]*hashCoeffs_[i];
   }
   if(M)
   {
      res = abs(res);
      if(hjump)
         *hjump = (res%(M-1))+1;
      res = res%M;
   }
   return res;
}


// Add a point to the LSH hash table using double hashing 
void FAMS::AddDataToHash(block HT[],int hs[],fams_point& pt,int where,int Bs,int M,int which,int which2,
                      int hjump)
{
   int nw=0;
   for(;;where = (where+hjump)%M)
   {
      nw++;
      if(nw > M)
      {
         bgLog("LSH hash table overflow exiting\n");
         exit(-1);
      }
      if(hs[where] == Bs)
         continue;
      HT[where][hs[where]].pt_ = &pt;
      HT[where][hs[where]].whichCut_ = which;
      HT[where][hs[where]].which2_ = which2;
      hs[where]++;
      break;
   }
}

// compute the pilot h_i's for the data points
void FAMS::ComputePilot(block *HT,int *hs, fams_partition *cuts, char* pilot_fn)
{
   const int win_j = 10,max_win=7000;
   int i,j;
   unsigned int nn;
   unsigned int wjd = (unsigned int)(win_j*d_);
   int num_l[1000];
   fams_res_cont res(n_);
   if(LoadBandwidths(pilot_fn)==0)
   {
      bgLog("compute bandwidths...");
      for(j=0; j<n_; j++)
      { 
         int numn=0;
         int numns[max_win/win_j];
         memset(numns,0,sizeof(numns));
         int nel;
         if(noLSH_)
         {
            nel = n_;
            for(i=0; i<nel; i++)
            {
               fams_pointp pt = &points_[i];
               nn = DistL1(points_[j],*pt) / wjd;
               if(nn <max_win/win_j)
                  numns[nn]++;	
            }
         }
         else
         {
            GetNearestNeighbours(points_[j],HT,hs,cuts,res,0,num_l);
            nel = res.nel_;      
            for(i=0; i<nel; i++)
            {
               fams_pointp pt = res.vec_[i];
               nn = DistL1(points_[j],*pt) / wjd;
               if(nn <max_win/win_j)
                  numns[nn]++;	
            }
         }
         for(nn=0; nn<max_win/win_j; nn++)
         {
            numn+=numns[nn];
            if(numn>k_)
            {
               break;
            }
         }	
         points_[j].window_=(nn+1)*wjd;
      }
      SaveBandwidths(pilot_fn);
   }
   else
      bgLog("load bandwidths...");
   for(j=0; j<n_; j++){
      points_[j].weightdp2_ = (float) pow(FAMS_FLOAT_SHIFT/points_[j].window_, (d_+2)*FAMS_ALPHA);
//      points_[j].weightdp2_ = (float) (1.0/pow(points_[j].window_, (d_+2)*FAMS_ALPHA));
   }
}

// compute real bandwiths for selected points
void FAMS::ComputeRealBandwidths(unsigned int h)
{
   const int win_j = 10,max_win=7000;
   int i,j;
   unsigned int nn;
   unsigned int wjd;
   wjd = (unsigned int) (win_j*d_);
   int who;
   if (h==0)
   {
      for(j=0; j<nsel_; j++)
      {
         who = psel_[j];
         int numn=0;
         int numns[max_win/win_j];
         memset(numns,0,sizeof(numns));
         for(i=0; i<n_; i++)
         {
            fams_pointp pt = &points_[i];
            nn = DistL1(points_[who],*pt) / wjd;
            if(nn <max_win/win_j)
               numns[nn]++;	
         }
         for(nn=0; nn<max_win/win_j; nn++)
         {
            numn+=numns[nn];
            if(numn>k_)
            {
               break;
            }
         }	
         points_[who].window_=(nn+1)*win_j;
      }
   } else
   {
      for(j=0; j<nsel_; j++)
      {
         who = psel_[j];
         points_[who].window_=h;
      }
   }
}

// compute the pilot h_i's for the data points
void FAMS::ComputeScores(block *HT,int *hs, fams_partition *cuts, float* scores)
{
   const int win_j = 10,max_win=7000;
   int i,j,who;
   unsigned int nn;
   unsigned int wjd = (unsigned int)(win_j*d_);
   int num_l[1000];
   memset(scores, 0, L_*sizeof(float));
   fams_res_cont res(n_);
   for(j=0; j<nsel_; j++)
   { 
      who = psel_[j];
      int numn;
      int nl=0;
      int numns[max_win/win_j];
      memset(numns,0,sizeof(numns));
      int nel;
      if(noLSH_)
      {
         nel = n_;
         num_l[L_]=n_+1;
         for(i=0; i<nel; i++)
         {
            fams_pointp pt = &points_[i];
            nn = DistL1(points_[who],*pt) / wjd;
            if(nn <max_win/win_j)
               numns[nn]++;	
            if(i == (num_l[nl]-1))
            {
               numn=0;
               for(nn=0; nn<max_win/win_j; nn++)
               {
                  if(numn>k_)
                     break;
               }
               for(; (num_l[nl]-1) == i; nl++)
                  scores[nl] += (float) (((nn+1.0)*win_j)/points_[who].window_);
            }
         }
      }
      else
      {
         GetNearestNeighbours(points_[who],HT,hs,cuts,res,0,num_l);
         nel = res.nel_;
         num_l[L_] = n_+1;
         for(i=0; i<nel; i++)
         {
            fams_pointp pt = res.vec_[i];
            nn = DistL1(points_[who],*pt) / wjd;
            if(nn <max_win/win_j)
               numns[nn]++;
            if(i == (num_l[nl]-1))
            {
               numn=0;
               for(nn=0; nn<max_win/win_j; nn++)
               {
                  numn+=numns[nn];
                  if(numn>k_)
                     break;
               }
               for(; (num_l[nl]-1) == i; nl++)
                  scores[nl] += (float) (((nn+1.0)*win_j)/points_[who].window_);
            }
         }
      }
      numn=0;
      for(nn=0; nn<max_win/win_j; nn++)
      {
         numn+=numns[nn];
         if(numn>k_)
            break;
      }	
   }
   for(j=0; j<L_; j++)
      scores[j]/=nsel_;

}

/*
Insert an mean-shift result into a second hash table so when another mean shift computationis 
performed about the same C_intersection region, the result can be retreived without further 
computation
*/

void FAMS::InsertIntoHash(block2* HT2,int *hs2,int where,int which,
                    unsigned short **solution,int M,int hjump)
{
   int nw=0;
   for(;;where = (where+hjump)%M)
   {
      nw++;
      if(nw==M)
      {
         bgLog("Second Hash Table Full\n");
         exit(-1);
      }
      if(hs2[where] == Bs2)
         continue;
      HT2[where][hs2[where]].dp_ = solution;
      HT2[where][hs2[where]++].whichCut_ = which;
      break;
   }
}

//perform an LSH query

void FAMS::GetNearestNeighbours(fams_point& who,block *HT, int *hs,
     fams_partition *cuts, fams_res_cont& res, int print,int num_l[])
{
   int i;
   for(i=0; i<L_; i++)
      EvalCutRes(who,cuts[i],t_cut_res_[i]);
   if(CompareCutRes(t_cut_res_,t_old_cut_res_)==0)
   {
      return;
   }
   memcpy(t_old_cut_res_,t_cut_res_,sizeof(t_old_cut_res_));
   res.clear();
   nnres1_++;
   for(i=0; i<L_; i++)
   {
      int hjump;
      int m = HashFunction(t_cut_res_[i],i,K_,M_,&hjump);
      int m2 = HashFunction(&t_cut_res_[i][1],i,K_-1);
      AddDataToRes(HT,hs,res,m,Bs,M_,i,nnres1_,m2,hjump);
      num_l[i] = res.nel_;
   }
}

// perform an LSH query using in addition the second hash table

unsigned short* FAMS::GetNearestNeighbours2H(unsigned short *who, block*HT, int *hs,
       fams_partition* cuts, fams_res_cont& res,
       unsigned short **solution, block2 *HT2, int *hs2)
{
   int i;
   for(i=0; i<L_; i++)
   {
      EvalCutRes(who,cuts[i],t_cut_res_[i]);
      t_m_[i] = HashFunction(t_cut_res_[i],i,K_,M_,&t_hjump_[i]);
      t_m2_[i] = HashFunction(&t_cut_res_[i][1],i,K_-1);
   }
   if(FAMS_DO_SPEEDUP)
   {
      int hjump2;
      int hf = HashFunction(t_m_,0,L_,M2_,&hjump2);
      int hf2 = HashFunction(&t_m_[L_/2-1],0,L_/2);
      unsigned short *old_sol =  FindInHash(HT2,hs2,hf,hf2,M2_,hjump2);
      if(old_sol!= NULL && old_sol != (unsigned short*)1)
         return old_sol;
      
      if(old_sol == NULL)
         InsertIntoHash(HT2,hs2,hf,hf2,solution,M2_,hjump2);
   }
   if(memcmp(t_m_,t_old_m_,sizeof(int)*L_)==0)
   {
      return NULL;
   }
   memcpy(t_old_m_,t_m_,sizeof(int)*L_);
   res.clear();
   nnres2_++;
   for(i=0; i<L_; i++)
      AddDataToRes(HT,hs,res,t_m_[i],Bs,M_,i,nnres2_,t_m2_[i],t_hjump_[i]);

   return NULL;
}

inline double SQ(double x) {return x*x;}

// perform an FAMS iteration

unsigned int FAMS::DoMeanShiftAdaptiveIteration(fams_res_cont& res,
              unsigned short *old, unsigned short *ret)
{
   double total_weight =0;
   int i,j;
   double dist;
   for(i=0; i<d_; i++) rr_[i]=0;
   int nel;
   if(noLSH_)
      nel = n_;
   else{
      nel = res.nel_;
   }
   fams_pointp ptp;
   unsigned int crtH;
   double hmdist=1e100;
   for(i=0; i<nel; i++)
   {
      ptp = noLSH_ ? &points_[i] : res.vec_[i];
      if(DistL1Data(old, *ptp, (*ptp).window_, dist))
      {
         double w;
         w = (*(ptp)).weightdp2_*SQ(1.0-(dist/(*ptp).window_));
         total_weight+=w;
         for(j=0; j<d_; j++)
            rr_[j] += (*(ptp)).data_[j]*w;
         if (dist<hmdist)
         {
            hmdist = dist;
            crtH=(*ptp).window_;
         }
      }
   }	
   if(total_weight==0)
   {
      return 0;
   }
   for(i=0; i<d_; i++)
      ret[i] = (unsigned short)(rr_[i]/total_weight);
   return crtH;
}

// perform a query on the second hash table

unsigned short* FAMS::FindInHash(block2* HT2,int *hs2,int where,int which,int M2,int hjump)
{
   int uu;
   int nw=0;
   for(;;where = (where+hjump)%M2)
   {
      nw++;
      if(nw > M2)
      {
         bgLog(" Hash Table2 full\n");
         exit(-1);
      }
      for(uu=0; uu<hs2[where]; uu++)
         if(HT2[where][uu].whichCut_ == which)
            return *(HT2[where][uu].dp_);
         if(hs2[where] < Bs2)
            break;
   }
   return NULL;
}


// Perform a query to one partition and retreive all the points in the cell 

void FAMS::AddDataToRes(block HT[],int hs[],fams_res_cont& res,int where,
        int Bs,int M,int which, unsigned short nnres, int which2, int hjump)
{
   int uu;
   for(;;where = (where+hjump)%M)
   {
      for(uu=0; uu<hs[where]; uu++)
      {
         if(HT[where][uu].whichCut_ == which && HT[where][uu].which2_ == which2)
         {
            if((*HT[where][uu].pt_).usedFlag_!=nnres)
            {
               res.push_back(HT[where][uu].pt_);
               (*HT[where][uu].pt_).usedFlag_=nnres;
            }
         }
      }
      if(hs[where] < Bs)
         break;
   }
}

// perform FAMS starting from a subset of the data points.
void FAMS::DoFAMS(block *HT,int *hs, fams_partition *cuts,
     block2* HT2, int *hs2)
{
   int jj;
   fams_res_cont res(n_);
   memset(HT2,0,sizeof(block2)*M2_);
   memset(hs2,0,sizeof(int)*M2_);
   unsigned short *oldMean;
   unsigned short *crtMean;
   oldMean = new unsigned short[d_];
   crtMean = new unsigned short[d_];
   fams_pointp currentpt;
   unsigned short *sol;
   //unsigned short *crtMode;
   int who;
   unsigned short **tMode;
   unsigned int newH;
   unsigned int *crtH;
   tMode = new unsigned short*[n_];
   bgLog(" Start MS iterations");
   int myPt = nsel_/10;
   for(jj=0; jj<nsel_; jj++)
   {
      if((jj%myPt)==0)
         bgLog(".");
      who = psel_[jj];
      currentpt = &points_[who];
      memcpy(crtMean, currentpt->data_, dataSize_);
      crtH = &hmodes_[jj];
      *crtH = currentpt->window_;
      tMode[jj]=(unsigned short*)1;
      int iter;
      for(iter=0; NotEq(oldMean,crtMean) && (iter<FAMS_MAXITER); iter++)
      {
         if(!noLSH_)
         {
            if((sol=GetNearestNeighbours2H(crtMean,HT,hs,cuts,
               res,&tMode[jj],
               HT2,hs2)))
            {
               if(sol == (unsigned short*)1)
               {
                  tMode[jj] = &modes_[jj*d_];
                  memcpy(tMode[jj], crtMean, dataSize_);
               }
               else
               {
                  tMode[jj] = &modes_[jj*d_];
                  memcpy(tMode[jj], sol, dataSize_);
                  break;
               }
            }
         }
         memcpy(oldMean, crtMean, dataSize_);
         if(!(newH=DoMeanShiftAdaptiveIteration(res,oldMean,crtMean)))
         {
            memcpy(crtMean, oldMean, dataSize_);
            break;
         }
         *crtH=newH;
      }
      if(tMode[jj]==(unsigned short*)1)
      {
         tMode[jj] = &modes_[jj*d_];
         memcpy(tMode[jj], crtMean, dataSize_);
      }
   }
   delete [] oldMean;
   delete [] crtMean;
   delete [] tMode;
   bgLog("done.\n");
}

void FAMS::SaveModes(char* fn)
{
   if (nsel_ < 1)
      return;
   bgLog("Save convergence points in %s ...", fn);
   FILE* fd;
   fd = fopen(fn, "wb");
   int i,j,idx;
   idx = 0;
   for (i=0; i<nsel_; i++)
   {
      for (j=0; j<d_; j++)
      {
         fprintf(fd, "%g ", modes_[idx++]*(maxVal_-minVal_)/65535.0 + minVal_);
      }
      fprintf(fd, "\n");
   }
   fclose(fd);
   bgLog("done\n");
}

void FAMS::SavePrunedModes(char* fn)
{
   if (npm_ < 1)
      return;
   bgLog("Save joined convergence points in %s ...", fn);
   FILE* fd;
   fd = fopen(fn, "wb");
   int i,j,idx;
   idx = 0;
   for (i=0; i<npm_; i++)
   {
      fprintf(fd, "%d  ", nprunedmodes_[i]);
      for (j=0; j<d_; j++)
      {
         fprintf(fd, "%g ", prunedmodes_[idx++]*(maxVal_-minVal_)/65535.0 + minVal_);
      }
      fprintf(fd, "\n");
   }
   fclose(fd);
   bgLog("done\n");
}
/*
int FAMS::PruneModes(int hprune, int npmin)
{
   bgLog(" Join Modes with h=%g and min pt=%d...", hprune*(maxVal_-minVal_)/65535.0, npmin);
   if (nsel_ < 1)
      return 1;
   hprune *= d_;

   int *mcount;
   float *cmodes, *ctmodes;
   unsigned short *pmodes;
   double cminDist, cdist;
   int iminDist, cref;
   mcount = new int[nsel_];
   cmodes = new float[d_*nsel_];

   int i, cd, cm, maxm;

   memset(mcount, 0, nsel_*sizeof(int));

   // set first mode
   for (cd = 0; cd<d_; cd++)
   {
      cmodes[cd] = modes_[cd];
   }
   mcount[0] = 1;
   maxm = 1;

   for (cm = 1; cm<nsel_; cm++)
   {
      bgLog("cm=%d, nm=%d\n",cm, maxm);

      pmodes = modes_+cm*d_;

      // compute closest mode
      cminDist = d_*1e7;
      iminDist = -1;
      for (cref = 0; cref<maxm; cref++)
      {
         cdist = 0;
         ctmodes = cmodes+cref*d_;
         for (cd=0; cd<d_; cd++)
            cdist += fabs(ctmodes[cd]/mcount[cref] - pmodes[cd]);
         if (cdist<cminDist)
         {
            cminDist = cdist;
            iminDist = cref;
         }
      }
      // join
      if (cminDist < hprune)
      {
         // aready in, just add
         for (cd=0; cd<d_; cd++)
         {
            cmodes[iminDist*d_+cd] += pmodes[cd];
         }
         mcount[iminDist] += 1;
      } else
      {
         // new mode, create
         for (cd=0; cd<d_; cd++)
         {
            cmodes[maxm*d_+cd] = pmodes[cd];
         }
         mcount[maxm] = 1;
         maxm += 1;
      }
   }

   // put the modes in the order of importance (count)
   int* stemp;
   int* istemp;
   stemp = new int[maxm];
   istemp = new int[maxm];
   for (i=0; i<maxm; i++)
   {
      stemp[i] = mcount[i];
      istemp[i] = i;
   }
   bgISort(stemp, maxm, istemp); // increasing

   // find number of relevant modes
   int nrel=1;
   for (i=maxm-2; i>=0; i--)
   {
      if (stemp[i]>=npmin)
         nrel++;
      else
         break;
   }

   CleanPrunedModes();
   prunedmodes_ = new unsigned short[d_*nrel];
   nprunedmodes_ = new int[nrel];
   unsigned short* cpm;
   npm_ = nrel;

   cpm = prunedmodes_;
   for (i=0; i<npm_; i++)
   {
      nprunedmodes_[i] = stemp[maxm-i-1];
      cm = istemp[maxm-i-1];
      for (cd=0; cd<d_; cd++)
      {
         *(cpm++) = (unsigned short) (cmodes[cm*d_+cd]/mcount[cm]);
      }
   }


   delete [] istemp;
   delete [] stemp;


   delete [] cmodes;
   delete [] mcount;

   bgLog("done\n");
   return 1;
}
*/

int FAMS::PruneModes(int hprune, int npmin)
{
   // compute jump
   int jm = ceil(((double) nsel_)/FAMS_PRUNE_MAXP);

   bgLog(" Join Modes with adaptive h/%d, min pt=%d, jump=%d\n", (int) pow(2,FAMS_PRUNE_HDIV), npmin, jm);
   bgLog("            pass 1");
   if (nsel_ < 1)
      return 1;
   hprune *= d_;

   int *mcount, *mcount2;
   float *cmodes, *ctmodes, *cmodes2;
   unsigned short *pmodes;
   double cminDist, cdist;
   int iminDist, cref;
   unsigned char *invalidm;
   invalidm = new unsigned char[nsel_];
   mcount = new int[nsel_];
   cmodes = new float[d_*nsel_];

   int i, cd, cm, maxm;

   memset(mcount, 0, nsel_*sizeof(int));
   memset(invalidm, 0, nsel_*sizeof(unsigned char));

   // set first mode
   for (cd = 0; cd<d_; cd++)
   {
      cmodes[cd] = modes_[cd];
   }
   mcount[0] = 1;
   maxm = 1;

   int myPt = FAMS_PRUNE_MAXP/10;
   for (cm = 1; cm<nsel_; cm+=jm)
   {
      if((cm%myPt)==0)
         bgLog(".");

      pmodes = modes_+cm*d_;

      //bgLog("cm=%d, nm=%d, kk=%d, %d %d\n",cm, maxm, hmodes_[cm], pmodes[0], pmodes[d_-1]);

      // compute closest mode
      cminDist = d_*1e7;
      iminDist = -1;
      for (cref = 0; cref<maxm; cref++)
      {
         if (invalidm[cref])
            continue;
         cdist = 0;
         ctmodes = cmodes+cref*d_;
         for (cd=0; cd<d_; cd++)
            cdist += fabs(ctmodes[cd]/mcount[cref] - pmodes[cd]);
         if (cdist<cminDist)
         {
            cminDist = cdist;
            iminDist = cref;
         }
      }
      // join
      hprune = hmodes_[cm] >> FAMS_PRUNE_HDIV;
      if (cminDist < hprune)
      {
         // aready in, just add
         for (cd=0; cd<d_; cd++)
         {
            cmodes[iminDist*d_+cd] += pmodes[cd];
         }
         mcount[iminDist] += 1;
      } else
      {
         // new mode, create
         for (cd=0; cd<d_; cd++)
         {
            cmodes[maxm*d_+cd] = pmodes[cd];
         }
         mcount[maxm] = 1;
         maxm += 1;
      }
      // check for valid modes
      if (maxm>2000)
      {
         for (i=0; i<maxm; i++)
         {
            if (mcount[i] < 3)
               invalidm[i] = 1;
         }
      }
   }

   bgLog("done\n");
   bgLog("            pass 2");
   
   // put the modes in the order of importance (count)
   int* stemp;
   int* istemp;
   stemp = new int[maxm];
   istemp = new int[maxm];
   for (i=0; i<maxm; i++)
   {
      stemp[i] = mcount[i];
      istemp[i] = i;
   }
   bgISort(stemp, maxm, istemp); // increasing

   // find number of relevant modes
   int nrel=1;
   for (i=maxm-2; i>=0; i--)
   {
      if (stemp[i]>=npmin)
         nrel++;
      else
         break;
   }
   if (nrel > FAMS_PRUNE_MAXM)
      nrel = FAMS_PRUNE_MAXM;

   // rearange only relevant modes
   mcount2 = new int[nrel];
   cmodes2 = new float[d_*nrel];

   for (i=0; i<nrel; i++)
   {
      cm = istemp[maxm-i-1]; // index
      mcount2[i] = mcount[cm];
      memcpy(cmodes2+i*d_, cmodes+cm*d_, d_*sizeof(float));
      //bgLog("1: %g %g %d\n",cmodes2[i*d_+0],cmodes2[i*d_+d_-1], mcount2[i]);
   }

   delete [] cmodes;
   memset(mcount, 0, nsel_*sizeof(int));
   mcount[0]=1;
   for (i=1; i<nsel_; i+=jm)
      mcount[i] = 1;

   maxm = nrel;

   myPt = nsel_/10;
   for (cm = 1; cm<nsel_; cm++)
   {
      if((cm%myPt)==0)
         bgLog(".");

      if (mcount[cm])
         continue;

      pmodes = modes_+cm*d_;

      // compute closest mode
      cminDist = d_*1e7;
      iminDist = -1;
      for (cref = 0; cref<maxm; cref++)
      {
         cdist = 0;
         ctmodes = cmodes2+cref*d_;
         for (cd=0; cd<d_; cd++)
            cdist += fabs(ctmodes[cd]/mcount2[cref] - pmodes[cd]);
         if (cdist<cminDist)
         {
            cminDist = cdist;
            iminDist = cref;
         }
      }
      // join
      hprune = hmodes_[cm] >> FAMS_PRUNE_HDIV;
      if (cminDist < hprune)
      {
         // aready in, just add
         for (cd=0; cd<d_; cd++)
         {
            cmodes2[iminDist*d_+cd] += pmodes[cd];
         }
         mcount2[iminDist] += 1;
      } else
      {
         // new mode, but discard in second pass
      }
   }

   // put the modes in the order of importance (count)
   for (i=0; i<maxm; i++)
   {
      stemp[i] = mcount2[i];
      istemp[i] = i;
   }
   bgISort(stemp, maxm, istemp); // increasing

   // find number of relevant modes
   nrel=1;
   for (i=maxm-2; i>=0; i--)
   {
      if (stemp[i]>=npmin)
         nrel++;
      else
         break;
   }

   CleanPrunedModes();
   prunedmodes_ = new unsigned short[d_*nrel];
   nprunedmodes_ = new int[nrel];
   unsigned short* cpm;
   npm_ = nrel;

   cpm = prunedmodes_;
   for (i=0; i<npm_; i++)
   {
      nprunedmodes_[i] = stemp[maxm-i-1];
      cm = istemp[maxm-i-1];
      for (cd=0; cd<d_; cd++)
      {
         *(cpm++) = (unsigned short) (cmodes2[cm*d_+cd]/mcount2[cm]);
      }
      //bgLog("2: %d %d\n",prunedmodes_[i*d_+0],prunedmodes_[i*d_+d_-1]);
   }


   delete [] istemp;
   delete [] stemp;

   delete [] cmodes2;
   delete [] mcount2;
   delete [] mcount;


   bgLog("done\n");
   return 1;
}

int FAMS::LoadBandwidths(char* fn)
{
   FILE* fd;
   fd = fopen(fn, "rb");
   if (fd == NULL)
      return 0;
   int n, i;
   fscanf(fd, "%d", &n);
   if (n!=n_)
   {
      fclose(fd);
      return 0;
   }
   float bw;
   float deltaVal = maxVal_-minVal_;
   for (i=0; i<n_; i++)
   {
      fscanf(fd, "%g", &bw);
      points_[i].window_ = (unsigned int) (65535.0*(bw)/deltaVal);
   }
   fclose(fd);
   return 1;
}

void FAMS::SaveBandwidths(char* fn)
{
   FILE* fd;
   fd = fopen(fn, "wb");
   if (fd==NULL)
      return;
   fprintf(fd, "%d\n", n_);
   float bw;
   float deltaVal = maxVal_-minVal_;
   int i;
   for (i=0; i<n_; i++)
   {
      bw = (float) (points_[i].window_*deltaVal/65535.0);
      fprintf(fd, "%g\n", bw);
   }
   fclose(fd);
}


// main function to find K and L
int FAMS::FindKL(int Kmin, int Kmax, int Kjump, int Lmax, int k, float width, float epsilon, int &K, int &L)
{
   bgLog("Find optimal K and L, K=%d:%d:%d, Lmax=%d, k=%d, Err=%.2g\n", Kmin, Kjump, Kmax, Lmax,
      k, epsilon);

   if (hasPoints_==0)
   {
      bgLog("Load points first\n");
      return 1;
   }

   int adaptive = 1;
   int hWidth=0;
   if (width>0)
   {
      adaptive = 0;
      hWidth = (int) (65535.0*(width)/(maxVal_-minVal_));
   }
   k_=k;
   epsilon += 1;

   // select points on which test is run
   SelectMsPoints(FAMS_FKL_NEL*100.0/n_, 0);

   // compute bandwidths for selected points
   ComputeRealBandwidths(hWidth);

   // start finding the correct l for each k
   float scores[FAMS_FKL_TIMES*FAMS_MAX_L];
   int Lcrt, Kcrt;
   
   int nBest;
   int LBest[FAMS_MAX_K];
   int KBest[FAMS_MAX_K];

   int ntimes, is;
   Lcrt = Lmax;
   bgLog(" find valid pairs");
   for (Kcrt = Kmax, nBest=0; Kcrt >= Kmin; Kcrt -= Kjump, nBest++)
   {
      // do iterations for crt K and L = 1...Lcrt
      for (ntimes=0; ntimes<FAMS_FKL_TIMES; ntimes++)
         DoFindKLIteration(Kcrt, Lcrt, &scores[ntimes*Lcrt]);
      
      // get correct for this k
      KBest[nBest]=Kcrt;
      LBest[nBest]=-1;
      for (is=0; (LBest[nBest]==-1) && (is<Lcrt); is++)
      {
         // find max on this column
         for (ntimes=1; ntimes<FAMS_FKL_TIMES; ntimes++)
         {
            if (scores[is]<scores[ntimes*Lcrt+is])
               scores[is] = scores[ntimes*Lcrt+is];
         }
         if (scores[is]<epsilon)
            LBest[nBest]=is+1;
      }
      bgLog(".");

      // update Lcrt to reduce running time!
      if (LBest[nBest]>0)
         Lcrt = LBest[nBest]+2;
   }
   bgLog("done\n");

   //start finding the pair with best running time
   double run_times[FAMS_FKL_TIMES];
   int iBest, i;
   double timeBest=-1;
   bgLog(" select best pair\n");
   for (i=0; i<nBest; i++)
   {
      if (LBest[i]<=0)
         continue;
      for (ntimes=0; ntimes<FAMS_FKL_TIMES; ntimes++)
         run_times[ntimes] = DoFindKLIteration(KBest[i], LBest[i], &scores[ntimes*Lcrt]);
      bgSort(run_times, FAMS_FKL_TIMES);
      if((timeBest==-1) || (timeBest>run_times[FAMS_FKL_TIMES/2]))
      {
         iBest = i;
         timeBest = run_times[FAMS_FKL_TIMES/2];
      }
      bgLog("  K=%d L=%d time: %g\n", KBest[i], LBest[i], run_times[FAMS_FKL_TIMES/2]);
   }
   K = KBest[iBest];
   L = LBest[iBest];

   bgLog("done\n");

   return 0;
}


double FAMS::DoFindKLIteration(int K,int L, float* scores)
{
   K_ = K;
   L_ = L;
   int i, j;

   // Allocate memory for the hash table
   M_ = GetPrime(3*n_*L_/(Bs));
   block *HT = new block[M_];
   int   *hs = new int[M_];
   InitHash(K_+L_);
   
   memset(hs,0,sizeof(int)*M_);

   // Build partitions
   fams_partition *cuts = new fams_partition[L_];
   for(i=0; i<20; i++)
      rand();
   int cut_res[FAMS_MAX_K];
   MakeCuts(cuts);

   //Insert data into partitions
   for(j=0; j<n_; j++)
   { 
      for(i=0; i<L_; i++)
      {
         EvalCutRes(points_[j],cuts[i],cut_res);
         int hjump;
         int m = HashFunction(cut_res,i,K_,M_,&hjump);
         int m2 = HashFunction(&cut_res[1],i,K_-1);
         AddDataToHash(HT,hs,points_[j],m,Bs,M_,i,m2,hjump);
      }
   }
   
   //Compute Scores
   timer_start();
   ComputeScores(HT, hs, cuts, scores);
   double run_time=timer_elapsed(0);

   // clean
   delete [] cuts;
   delete [] hs;
   delete [] HT;

   return run_time;
}

// main function to run FAMS
int FAMS::RunFAMS(int K, int L, int k, double percent, int jump, float width, char* pilot_fn)
{

   bgLog("Running FAMS with K=%d L=%d\n", K, L);
   if (hasPoints_==0)
   {
      bgLog("Load points first\n");
      return 1;
   }
   int i, j;

   int adaptive = 1;
   int hWidth;
   if (width>0)
   {
      adaptive = 0;
      hWidth = (int) (65535.0*(width)/(maxVal_-minVal_));
   }

   K_=K; L_=L; k_=k;
   SelectMsPoints(percent, jump);

   // Allocate memory for the hash table
   M_ = GetPrime(3*n_*L_/(Bs));
   M2_ = GetPrime((int)(nsel_*20*3/Bs2));
   block *HT = new block[M_];
   int   *hs = new int[M_];
   block2 *HT2 = new block2[M2_];
   int   * hs2 = new int[M2_];
   
   memset(hs,0,sizeof(int)*M_);

   // Build partitions
   fams_partition *cuts = new fams_partition[L];
   for(i=0; i<20; i++)
      rand();
   int cut_res[FAMS_MAX_K];
   MakeCuts(cuts);

   InitHash(K_+L_);

   //Insert data into partitions
   for(j=0; j<n_; j++)
   { 
      for(i=0; i<L_; i++)
      {
         EvalCutRes(points_[j],cuts[i],cut_res);
         int hjump;
         int m = HashFunction(cut_res,i,K_,M_,&hjump);
         int m2 = HashFunction(&cut_res[1],i,K_-1);
         AddDataToHash(HT,hs,points_[j],m,Bs,M_,i,m2,hjump);
      }
   }	

   //Compute pilot if necessary
   bgLog(" Run pilot ");
   if(adaptive)
   {
      bgLog("adaptive...");
      ComputePilot(HT,hs,cuts,pilot_fn);
   }
   else
   {
      bgLog("fixed bandwith...");
      unsigned int hwd = (unsigned int)(hWidth*d_);
      for(i=0; i<n_; i++)
      {
         points_[i].window_=hwd;
         points_[i].weightdp2_=1;
      }
   }
   bgLog("done.\n");

   DoFAMS(HT, hs, cuts, HT2, hs2);

   // join modes
   PruneModes(FAMS_PRUNE_WINDOW, FAMS_PRUNE_MINN);

   // clean
   delete [] cuts;
   delete [] hs2;
   delete [] HT2;
   delete [] hs;
   delete [] HT;

   return 0;
}

void usage()
{
   bgLog("usage: K L k file_name input_dir [ [-j jump | -p percent] | -h width \n   -f epsilon Kmin Kjump]"); 
}


int main(int argc,char** argv)
{
   if(argc < 6)
   {
      usage();
      exit(1);
   }
   int no_lsh, find_kl;
   int K, L, k_neigh;
   char *data_file_name, *input_path;
   char fdata_file_name[200];
   K = atoi(argv[1]);
   L = atoi(argv[2]);
   no_lsh = (K<=0) || (L <= 0);
   find_kl = 0;

   float epsilon;
   int Kmin, Kjump, Kmax;
   int Lmax;

   float width = -1;

   k_neigh = atoi(argv[3]);
   data_file_name = argv[4];
   input_path = argv[5];
   sprintf(fdata_file_name,"%s%s.txt",input_path,data_file_name);
   //printf("Usage: FastAdaptiveMeanShift_2D <K> <L> <k> <input text file with per-entry multivariate ushort-range values> <directory of input- and temporary files>");
   //printf("K - 'observation size'");
   //printf("L - dependable size");
   //printf("k - binning size for LSH tree");
   int jump=1;
   double percent=0.0;
   int i;
   if (argc > 6)
   {
      for (i=6; i<argc; i++)
      {
         if (argv[i][0] != '-')
         {
            bgLog("Error in param %s\n", argv[i]);
            usage();
            exit(1);
         }
         switch(argv[i][1])
         {
         case 'j': 
            i++;
            jump = atoi(argv[i]);
            if (jump<1) jump=1;
            break;
         case 'p':
            i++;
            percent = atof(argv[i]);
            if ((percent<0) || (percent>1)) percent = 0;
            break;
         case 'h':
            i++;
            width = (float) atof(argv[i]);
            break;
         case 'f':
            i++;
            epsilon = (float) atof(argv[i++]);
            Kmin = atoi(argv[i++]);
            Kjump = atoi(argv[i]);
            Lmax = L; Kmax = K;
            find_kl=1;
            break;
         default:
            bgLog("Error in param %s\n", argv[i]);
            usage();
            exit(1);
            break;
         }
      }
   }


   // load points
   FAMS cfams(no_lsh);
   if (cfams.LoadPoints(fdata_file_name))
      return 1;

   // find K L (if necessary)
   if (find_kl)
   {
      cfams.FindKL(Kmin, Kmax, Kjump, Lmax, k_neigh, width, epsilon, K, L);
      bgLog("Found K = %d L = %d (write them down)\n", K, L);
   //   int ch=' ';
   // DISABLING STRAIGHT COMPUTATION TO AVOID INTERACTIVE MODE
   //   do{
   //      bgLog("Do you want to run FAMS with this (K=%d,L=%d) pair? (y/n)",K,L);
   //      ch = getchar();
   //      if ((ch == 'n') || (ch == 'N'))
   //         return 0;
   //   } while((ch != 'y') && (ch != 'Y'));
      return 0;
   }
   sprintf(fdata_file_name, "%spilot_%d_%s.txt", input_path, k_neigh, data_file_name);
   cfams.RunFAMS(K, L, k_neigh, percent, jump, width, fdata_file_name);

   // save the data
   sprintf(fdata_file_name,"%sout_%s.txt",input_path,data_file_name);
   cfams.SaveModes(fdata_file_name);

   // save pruned modes modes
   sprintf(fdata_file_name,"%smodes_%s.txt",input_path, data_file_name);
   cfams.SavePrunedModes(fdata_file_name);


   return 0;
}

