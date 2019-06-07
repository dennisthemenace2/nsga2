
#include <iostream>
#include <vector>
#include <list>

#include <math.h>

#include <random>
    

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <string>

#include <stdio.h>
#include <bitset>


#include <pthread.h>


#define INF 1.0e14

using namespace std;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)
  
  // For more on using Rcpp click the Help button on the editor toolbar

#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

//nsga2(simulateBasicModel, 
          //    GASettingsList$parameterNr, 
          //    GASettingsList$objectiveNr,
        //      generations=GASettingsList$generationNr,
        //      popsize=GASettingsList$popSize, # was 12, 12 
        //      cprob=GASettingsList$cProb,  
      //        cdist=GASettingsList$cDist,   
    //          mprob=GASettingsList$mProb,  
    //          mdist=GASettingsList$mDist,   
    //          lower.bounds=rep(0, GASettingsList$parameterNr),
    //          upper.bounds=rep(1, GASettingsList$parameterNr))
              
              

///function(generations, popsize, cros,mutation,params,fn,objectives, boundmin=0, boundmax=1  ){
 
struct nsgaConfig{
   unsigned int nparams;
   unsigned int nobj;
       
   double lowerbound; 
   double upperbound;
   double crossrate;
   double eta_c;
   double mutationrate;
   double eta_m;
       
};
     
struct individual{
    int rank;
    double crowd_dist;

    std::vector<double> err;
    std::vector<double> params;
    
    int domcnt;
    std::list<individual*> domset;
    
};


class cmp {

public:
    cmp(unsigned p) : idx(p) {}
    unsigned idx;

    bool operator()(individual *p, individual *q) {
        return p->err[idx]< q->err[idx];
    }
};

bool sortBydist( individual* p,  individual* q) { return p->crowd_dist > q->crowd_dist; }


void initPopulation(std::vector<individual> &population, std::vector<individual*> &parents,  std::vector<individual*> &replace, nsgaConfig *cfg){
 // unsigned int BITS = sizeof(unsigned int); // could be customized
  
//  unsigned int needed = ceil(  (nparams * BITS) / sizeof(unsigned int) );
  unsigned int needed = cfg->nparams;
  unsigned psize  =parents.size();
  
  for(unsigned i=0;i<population.size();++i ){  
    // init population.
    population[i].params.resize(needed) ;
    population[i].err.resize(cfg->nobj);
     
    for(unsigned k=0;k<needed;++k ){
      double delta = cfg->upperbound - cfg->lowerbound;
      population[i].params[k] = cfg->lowerbound + delta*(rand() / double(RAND_MAX) );
    }

    if(i<psize){
      parents[i]= &population[i]; // init parnts
    }else{
      replace[i-psize]= &population[i]; // init replace 
    }
    
  }
  
}


///THREAD



void reportCurrentPopErr(std::vector<individual*> &parents,unsigned nobj,  Rcpp::Function func){

    Rcpp::NumericMatrix currerr(parents.size(), nobj);
    
    for(unsigned i=0;i<parents.size();++i ){      
      for(unsigned a=0;a<nobj;++a){
        currerr(i,a) = parents[i]->err[a];
      }    
    }
    func( currerr );

}

void evaluateVec(std::vector<individual*> &population, Rcpp::Function func,nsgaConfig *cfg){
  
   Rcpp::NumericMatrix numparams(population.size(),cfg->nparams);

  //prepare matrix
  for(unsigned i=0;i<population.size();++i ){      
    for(unsigned a=0;a<cfg->nparams;++a){
      numparams(i,a) = population[i]->params[a];
    }    
  }
  
  Rcpp::NumericMatrix reterr = func( numparams );
  
// check return real quick
  if(reterr.nrow() != population.size() ){
    cout<<"matrix in wrong format"<<endl;
    cout<<"especting:"<< population.size()<<"x" <<cfg->nobj <<endl;
    cout<<"got:"<< reterr.nrow()<<"x" <<reterr.ncol() <<endl;
    Rcpp::stop("Wrong matrix dimensions!\n"); 
  }
  
  //write back results
  for(unsigned i=0;i<population.size();++i ){      
    for(unsigned a=0;a<reterr.ncol();++a){
        population[i]->err[a]= reterr(i,a);
    }    
  }
   
}

void evaluate(std::vector<individual*> &population, Rcpp::Function func){
  
  for(unsigned i=0;i<population.size();++i ){ 
     Rcpp::NumericVector err =  func( population[i]->params );
     for( unsigned k=0;k<err.size();++k ){ // not good but well
       population[i]->err[k]= err[k]; //copy
     }
  }
  
}


bool dominates(individual& p, individual& q ){
  
  for(unsigned x=0; x< p.err.size();++x){
    if(p.err[x]>q.err[x] ){
      return false;
    }
  }
  
  return true;
}

void fastNonDominatedSort(std::vector<individual> &population,unsigned from, unsigned to,std::list<std::list<individual*> > &ret ){
  
  ret.clear();
  
  std::list<individual*> firstfront;
  for( unsigned x=from; x<to;++x){
    population[x].domset.clear();
    population[x].domcnt=0;
    
    for(unsigned y=from; y<to;++y){
      if( dominates(population[x],population[y]) ){
        population[x].domset.push_back(&population[y]);
        //cout<<"pushback donination"<<"\n";
      }else if(dominates(population[y],population[x]) ){
        ++population[x].domcnt;
      } 
    }
    if(population[x].domcnt==0){
      population[x].rank = 1 ;
   //  cout<<"add first front\n";
      firstfront.push_back(&population[x]);
    }
    
  }
  
// cout<<"first front size:"<<firstfront.size()<<"\n";
  ret.push_back(firstfront);
  
  std::list<individual*> *nextfront = &firstfront;
 
  //std::list<individual*> nextfront = firstfront;
 //
 // std::list<individual*> *domfront;
    
  std::list<individual*> addfront;
 // std::list<individual*> nfront;
    
  std::list<individual*>::iterator it;
  std::list<individual*>::iterator it2;
    
  unsigned curr=2;
 // cout<<"nextfront size;"<<nextfront->size()<<"\n";
  while(nextfront->size()  >0 ){

    addfront.clear();

    for (it=nextfront->begin(); it!=nextfront->end(); ++it){
   //   cout<<"(*it)->domset.size()"<<(*it)->domset.size()<<"\n";
        for(it2 = (*it)->domset.begin() ; it2 != (*it)->domset.end();++it2 ){
          --((*it2)->domcnt);
          if((*it2)->domcnt == 0){
            (*it2)->rank = curr;
            addfront.push_back(*it2);
           // cout<<"add front "<< curr<<"\n";
          //  nextfront.push_back(*it2);
          }
          
        }
    }//end for
    if(addfront.empty() ){ // because back on empty list is undefined
      break;
    }
    ++curr;
    ret.push_back(addfront);
    //nextfront = addfront;
  //   cout<<"pushback next front size:"<<addfront.size()<<"\n";
     
    nextfront = &ret.back(); 
  
  }//end while


}


void crossover(individual *p ,individual *q,individual* replace,individual* replace2, nsgaConfig *cfg){
  
  if( (rand() / double(RAND_MAX) ) < cfg->crossrate){
    replace->params = p->params;
    return;
  }
  double y1,y2;
  double betaq;
  
  for(unsigned i =0; i< p->params.size() ;++i){
    if( (rand()/double(RAND_MAX)) <= 0.5 ){
        if ( p->params[i] < q->params[i]) {
          y1 = p->params[i];
          y2 = q->params[i];
        } else {
          y1 = q->params[i];
          y2 = p->params[i];
        }
        double yl =  cfg->lowerbound;  // for individual bounds
        double yu = cfg->upperbound;
        double randv = (rand()/double(RAND_MAX));
        double beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
        double alpha = 2.0 - pow(beta,-(cfg->eta_c+1.0));
        if (randv <= (1.0/alpha)){
          betaq = pow ((randv*alpha),(1.0/(cfg->eta_c+1.0)));
        }else{
          betaq = pow ((1.0/(2.0 - randv*alpha)),(1.0/(cfg->eta_c+1.0)));
        }
        double c1 = 0.5*((y1+y2)-betaq*(y2-y1));
        beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
        alpha = 2.0 - pow(beta,-(cfg->eta_c+1.0));
        if (randv <= (1.0/alpha)){
          betaq = pow ((randv*alpha),(1.0/(cfg->eta_c+1.0)));
        }else{
          betaq = pow ((1.0/(2.0 - randv*alpha)),(1.0/(cfg->eta_c+1.0)));
        }
        double c2 = 0.5*((y1+y2)+betaq*(y2-y1));
        /* Enforce constraints: */
        if (c1 < yl) c1=yl;
        if (c2 < yl) c2=yl;
        if (c1 > yu) c1=yu;
        if (c2 > yu) c2=yu;
        if ( (rand() / double(RAND_MAX) ) <= 0.5) {
          replace->params[i] = c2;
          replace2->params[i] = c1;
        } else  {
          replace->params[i] = c1;
          replace2->params[i] = c2;
        }
        
    }else{// donot change
      replace->params[i] =  p->params[i];
      replace2->params[i] =  q->params[i];
    }
    

  }
    
}
void mutation(individual *p , nsgaConfig *cfg){
  
  double val, xy,deltaq;
  
  for(unsigned i =0; i< p->params.size() ;++i){
    
    if( (rand() / double(RAND_MAX)) < cfg->mutationrate){
      double y = p->params[i];
      
      double yl = cfg->lowerbound;
      double yu = cfg->upperbound;
      double delta1 = (y-yl)/(yu-yl);
      double delta2 = (yu-y)/(yu-yl);
      double randv =  (rand() / double(RAND_MAX)) ;
      
      double mut_pow = 1.0/(cfg->eta_m+1.0);
      if(randv <= 0.5){
        xy = 1.0-delta1;
        val = 2.0*randv+(1.0-2.0*randv)*(pow(xy,(cfg->eta_m+1.0)));
        deltaq =  pow(val,mut_pow) - 1.0;
      }else{
        xy = 1.0-delta2;
        val = 2.0*(1.0-randv)+2.0*(randv-0.5)*(pow(xy,(cfg->eta_m+1.0)));
        deltaq = 1.0 - (pow(val,mut_pow));
      }
      y = y + deltaq*(yu-yl);
      if (y < yl){ 
        y = yl;
      }
      if (y > yu){
        y = yu;
      }
      p->params[i] = y;
        
    }
      //else skip
  }
  
}
void reproduce(std::vector<individual*> &selected, std::vector<individual*> &replace, nsgaConfig *cfg){

  for(unsigned i =0;i<selected.size();i+=2 ){
    
  //crossover( selected[uni(rng)], selected[uni(rng)] ,replace[i],replace[i+1] , cfg);
    crossover( selected[i], selected[i+1] ,replace[i],replace[i+1] , cfg);
    mutation(replace[i],cfg);/// use etac for test !!!
  }

}

void calculate_crowding_distance(std::list<individual*>  &front,unsigned objectives){
  unsigned N = front.size();
  //cout<<"calculate_crowding_distance fronts.size()"<<front.size()<<"\n";
  
  if(N <= 2){
    front.front()->crowd_dist = INF;
    if(N == 1){
      return;
    }
    
    front.back()->crowd_dist = INF;
    return;
  }
  
  std::vector<individual*> frontvec(N);
  std::list<individual*>::iterator it;
  
  unsigned idx =0; // copy pointers to vector.
  for (it=front.begin(); it!=front.end(); ++it){
    frontvec[idx++] = (*it);
    (*it)->crowd_dist=0;
  }
  cmp compr( 0);
  for(unsigned i =0;i<objectives;++i){
    compr.idx=i;
    std::sort(frontvec.begin(), frontvec.end(), compr );
  //  cout<<"now sorted\n";
    
  /*  for(unsigned w = 0;w<frontvec.size();++w ){
     cout<<"sorted vec:"<< frontvec[w]->err[i]<<"\n";
    }*/
    
    individual* tmp  = frontvec.back();
    tmp->crowd_dist = INF;
    double reg = tmp->err[i];     
    tmp = frontvec.front();
    tmp->crowd_dist = INF;
    reg -= tmp->err[i];

    for(unsigned k=1;k<(N-1);++k){
      if (frontvec[k]->crowd_dist != INF){    
        frontvec[k]->crowd_dist +=  (frontvec[k+1]->err[i]-frontvec[k-1]->err[i]) / reg;
      }
      //else{
    //   cout<<"ok was inf alrady\n";
    //  }
    }
    
  }
  
  //thi sis debug only
//   for(unsigned w = 0;w<frontvec.size();++w ){
//     cout<<"cdist of vec:"<< frontvec[w]->crowd_dist<<"\n";
//    }
 
}

void select_parents(std::list<std::list<individual*> > &fronts, unsigned objectives,std::vector<individual*> &parents,std::vector<individual*> &replace){
  
  unsigned inv = 0;
  unsigned Npop = parents.size();

   std::list<std::list<individual*> >::iterator it;
   std::list<individual*>::iterator it2;
  
  //cout<<"fronts to select from"<<fronts.size() << " \n";
  
   for (it=fronts.begin(); it!=fronts.end(); ++it){
   //  cout<<"for all fronts\n";
     
     calculate_crowding_distance(*it,objectives);
  //  cout<<"crowdin dist is assigned\n";
     
     if( (*it).size() +inv > Npop ) {
       //select rest..
    //   cout<<"prepare rest\n";
      
       std::vector<individual*> indvVec( (*it).size() );
       unsigned cnt=0;
       for (it2=(*it).begin(); it2!= (*it).end(); ++it2){
         indvVec[cnt++]= (*it2);
       }
       std::sort(indvVec.begin(), indvVec.end(), sortBydist );
       
     //  cout<<"vector is sorted\n";
      // for(unsigned w=0;w<indvVec.size();++w){
    //     cout<<"sorted vec by crowd dist"<<indvVec[w]->crowd_dist<<"\n";
    //   }
       
       // vector now sorted
       unsigned repcnt = 0;
       for( unsigned k=0;k<indvVec.size();++k ){
         if(inv < Npop){
           parents[inv++]= indvVec[k];
         }else{
           replace[repcnt++] = indvVec[k];
         }
       }
    //   cout<<"append rest"<<repcnt<<"\n";
       
       for (++it; it!=fronts.end(); ++it){ // set replace list
          for (it2=(*it).begin(); it2!= (*it).end(); ++it2){
            replace[repcnt++]= (*it2);
        //    cout<<"append rest cnt"<<repcnt<<"\n";
             
          }
       }
       
       break;
     }
     
     //cout<<"append whole front"<<inv<<"\n";
        
     for (it2=(*it).begin(); it2!= (*it).end(); ++it2){
        parents[inv++]= (*it2);
    //    cout<<"append surviver"<<inv<<" rank:"<< (*it2)->rank<<"\n";
        //   cout<<"look what parent is now:"<<parents[inv-1]->rank<<"\n";
     }

   }
   
          

}



void tournament(std::vector<individual*> &parents, std::vector<individual*> &selected,std::uniform_int_distribution<int> &uni,std::mt19937 &rng){
  

  //tournement selection sould be vectorized somhow
  for( unsigned i=0;i<parents.size();++i ){
    unsigned a = uni(rng);
    unsigned b = uni(rng);
 //   unsigned a = i;
 //     unsigned b = (i+1)*(i< parents.size()-1) ;
    
    if(parents[a]->rank < parents[b]->rank){
      selected[i] = parents[a];
      //    cout<<"selected:"<<a<<"\n";
    }else if(parents[a]->rank > parents[b]->rank){
      selected[i] = parents[b];
      //     cout<<"selected:"<<b<<"\n";
      
    }else{
      if(parents[a]->crowd_dist > parents[b]->crowd_dist){
        selected[i] = parents[a];
        //    cout<<"selected:"<<a<<"\n";
        
      }else{
        selected[i] = parents[b];
        //    cout<<"selected:"<<b<<"\n";
        
      }
    }
  }
  
}

 
// [[Rcpp::export]]
Rcpp::List my_nsga2(Rcpp::Function func,Rcpp::NumericVector params,Rcpp::NumericVector objectives ,
                    Rcpp::NumericVector generations,Rcpp::NumericVector popsize,Rcpp::NumericVector cros,
                    Rcpp::NumericVector crosIdx,
                    Rcpp::NumericVector mutation,Rcpp::NumericVector mutationIdx, Rcpp::NumericVector lower, Rcpp::NumericVector upper){
    
    // shoudl check input
    if( params.size() != 1 ||popsize.size() != 1 ){
      cout<<"No value for params set\n";
    }
    unsigned int psize = popsize[0];
    unsigned int ngenerations =generations[0];
    
   
    //chk popsize real quick
    if(psize % 2 != 0){
      cout<<"Population size must be multible of 2\n\n";
  //   Rcpp::stop("Population size must be multible of 2\n"); 
      return Rcpp::List::create();
    }
//  unsigned int nparams = params[0];
//    unsigned int nobj = objectives[0];
 //   double lowerbound = lower[0]; // fix the bound
  //  double upperbound = upper[0];
  //  double crossrate = cros[0];
//    double mutationrate = mutation[0];
    
    std::list<std::list<individual*> > fronts ;
       
    std::vector<individual> population(psize*2) ;
    
    nsgaConfig cfg;
    cfg.nparams = params[0];
    cfg.nobj = objectives[0];
    
    cfg.crossrate =cros[0];
    cfg.eta_c=crosIdx[0];
    
    cfg.lowerbound=lower[0]; 
    cfg.upperbound=upper[0];
    
    cfg.mutationrate=mutation[0];
    cfg.eta_m=mutationIdx[0];
    
    
    
  // set up config for easy access
  

    
 //   cout<<"start here\n";
    
    std::vector<individual*> parents(psize); // these to I need
    std::vector<individual*> replace(psize);    
    
    initPopulation(population,parents,replace,&cfg); // allocates memory
  //  cout<<"population inited\n";

    
     //evaluate(parents, func);
    evaluateVec(parents, func,&cfg);
  //       cout<<"evaluate!\n";

    fastNonDominatedSort(population,0, psize,fronts);
 //   cout<<"fast non dom sort!\n";
    // select
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,(psize-1) ) ; // guaranteed unbiased
 
    std::vector<individual*> selected(psize);


    tournament(parents, selected,uni,rng);          
    reproduce(selected,replace, &cfg);
    //evaluate(replace, func);
    evaluateVec(replace, func,&cfg);

    unsigned P2 =  psize*2;
    for(unsigned i =0;i<ngenerations;++i){
      //population is already union of both...
       fastNonDominatedSort(population,0,P2,fronts);
       
  //     cout<<"frons estimated :"<<fronts.size()<<"\n";
       select_parents(fronts,cfg.nobj,parents ,replace);
       
    //cout<<"select parents\n";
  
       tournament(parents, selected,uni,rng);                  
       reproduce(selected,replace, &cfg);
       
    //   evaluate(replace,func);     
       evaluateVec(replace,func, &cfg);
       
    }
    fastNonDominatedSort(population,0,P2,fronts);
    select_parents(fronts,cfg.nobj,parents ,replace); 
    // parents has now the results.
//    for(unsigned i =0;i<parents.size();++i){
//      cout<<"lokk rank:"<<parents[i]->rank<<"\n";
//    }
    
    Rcpp::NumericMatrix numerr(psize,cfg.nobj );
    Rcpp::NumericMatrix numparams(psize,cfg.nparams );
    
// std::vector<double> decparams(nparams);
    std::vector<int> frank(psize);
    
    
    //prepare results.
    for(unsigned k=0;k<psize;++k){
//      decode( parents[k]->params , lowerbound, upperbound,decparams);
      frank[k] = parents[k]->rank;
   //   cout<<"rank "<<parents[k]->rank<<"\n";
      
      for(unsigned a=0;a<cfg.nparams;++a){
        numparams(k,a) =  parents[k]->params[a];
      }    
      for(unsigned a=0;a<cfg.nobj;++a){
        numerr(k,a) = parents[k]->err[a];
      }    
      
    }
    
    //deletePolpulation(population); // free memory
    return Rcpp::List::create(Rcpp::Named("error") = numerr,
                              Rcpp::Named("parameters") = numparams,
                              Rcpp::Named("rank") = frank
                               
                              );
}

