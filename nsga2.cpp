
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

void decode( std::vector<unsigned> &params,double boundmin,  double boundmax,std::vector<double> &result){
 
  unsigned nparams = params.size();
  double t = pow(2,sizeof(unsigned)*8); //(1<<(sizeof(unsigned)*8)-1);     
  double w = (boundmax- boundmin);
  
  for(unsigned i=0;i<nparams;++i){
    result[i] = boundmin +  params[i] * (w / (double)(t-1)) ;
  }

}
    
     
struct individual{
    int rank;
    double crowd_dist;


    std::vector<double> err;
    std::vector<unsigned> params;
    
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


void initPopulation(unsigned int nparams, unsigned objectives, std::vector<individual> &population, std::vector<individual*> &parents,  std::vector<individual*> &replace){
 // unsigned int BITS = sizeof(unsigned int); // could be customized
  
//  unsigned int needed = ceil(  (nparams * BITS) / sizeof(unsigned int) );
  unsigned int needed = nparams;
  unsigned psize  =parents.size();
  
  for(unsigned i=0;i<population.size();++i ){  
    // init population.
    population[i].params.resize(needed) ;
    population[i].err.resize(objectives);
     
    for(unsigned k=0;k<needed;++k ){
      population[i].params[k]= (rand()<<16) | (rand()& 0xffff) ;// not evnough bits...
      
      //std::bitset<32>      x(rand()<<16);
    //  std::cout << x << std::endl;
      
    }
  /*  std::vector<double> decres(nparams);
    decode(population[i].params, -10,10,decres);
    
    for(unsigned w = 0;w<nparams;++w){
      cout<<"decoded:"<<decres[w]<<"\n";
    }*/
    
    if(i<psize){
      parents[i]= &population[i]; // init parnts
    }else{
      replace[i-psize]= &population[i]; // init replace 
    }
    
  }
  
}


///THREAD


struct evalThread_data{
   unsigned  from;
   unsigned  to;
   unsigned  nparams;
   double lower;
   double uper;
   
//   std::vector<double> decodeprams; // pass nparams
   std::vector<individual*> *population;     
   Rcpp::Function *func;
   
};


void *evalThread(void *arg){
   
   struct evalThread_data *data;
   data = (struct evalThread_data *) arg;
   
   cout<<"welcome to thread: create decode vector"<< data->nparams<<"\n";
   std::vector<double> decodeprams(data->nparams);
   
   for(unsigned i=data->from;i<data->to;++i ){ 
     
     decode( (*data->population)[i]->params , data->lower,data->uper,decodeprams );
     cout<<"decode done"<<i<<"\n";
     
     Rcpp::NumericVector err =  (*data->func)( decodeprams );
  
     for( unsigned k=0;k<data->nparams;++k ){
       (*data->population)[i]->err[k]= err[k]; //copy
     //  cout<<"k"<<k<<" "<< err[k]<<"\n";
     }
   
  }
  
    cout<<"thread done\n";

   pthread_exit( (void *)NULL);
}

void evaluateMT(std::vector<individual*> *parents ,std::vector<evalThread_data> &tdata,std::vector<pthread_t> &thread,pthread_attr_t *attr){
  // do calculations.
      cout<<"prepare"<<"\n";
      
  // start threads
       for(unsigned t=0; t<tdata.size(); ++t)  {   
           cout<<"start thread"<<"\n";
               
           tdata[t].population = parents;// set target
           pthread_create(&thread[t],attr, evalThread, (void *)&tdata[t]);
       }

      // wait for all threads
       for(unsigned t=0; t<tdata.size(); ++t)  {
           cout<<"wait for thread end"<<"\n";
          pthread_join(thread[t], NULL);
       }
       cout<<"evaluate done"<<"\n";
}


void reportCurrentPopErr(std::vector<individual*> &parents,unsigned nobj,  Rcpp::Function func){

    Rcpp::NumericMatrix currerr(parents.size(), nobj);
    
    for(unsigned i=0;i<parents.size();++i ){      
      for(unsigned a=0;a<nobj;++a){
        currerr(i,a) = parents[i]->err[a];
      }    
    }
    func( currerr );

}
void evaluateVec(std::vector<individual*> &population,unsigned nparams,  Rcpp::Function func,double lower,double upper,unsigned nobj){
  
   std::vector<double> decodeprams( nparams ); // pass nparams

   Rcpp::NumericMatrix numparams(population.size(),nparams );

  //prepare matrix
  for(unsigned i=0;i<population.size();++i ){      
    decode( population[i]->params , lower,upper,decodeprams );
    for(unsigned a=0;a<nparams;++a){
        numparams(i,a) = decodeprams[a];
    }    
  }
  
  Rcpp::NumericMatrix reterr = func( numparams );
  
// check return real quick
  if(reterr.nrow() != population.size() ){
    cout<<"matrix in wrong format"<<endl;
    cout<<"especting:"<< population.size()<<"x" <<nobj <<endl;
    cout<<"got:"<< reterr.nrow()<<"x" <<reterr.ncol() <<endl;
    Rcpp::stop("Wrong matrix dimensions!\n"); 
  }
  //write back results
  for(unsigned i=0;i<population.size();++i ){      
    decode( population[i]->params , lower,upper,decodeprams ); // can i remove this line
    for(unsigned a=0;a<reterr.ncol();++a){
        population[i]->err[a]= reterr(i,a);
     //  cout<<"re err"<<population[i]->err[a];
    }    
   //  cout<<"\n";
        
  }
   
}

void evaluate(std::vector<individual*> &population,unsigned nparams,  Rcpp::Function func,double lower,double upper){
  
  std::vector<double> decodeprams( nparams ); // pass nparams

  for(unsigned i=0;i<population.size();++i ){ 
     
     decode( population[i]->params , lower,upper,decodeprams );
     
     Rcpp::NumericVector err =  func( decodeprams );
  
     for( unsigned k=0;k<err.size();++k ){ // not good but well
       population[i]->err[k]= err[k]; //copy
     //  cout<<"k"<<k<<" "<< err[k]<<"\n";
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


void crossover(individual *p ,individual *q,individual* replace, double crossrate){
  
  if( (rand() / double(RAND_MAX) ) < crossrate){
    replace->params = p->params;
    return;
  }
  
  for(unsigned i =0; i< p->params.size() ;++i){
    int k =  (rand()<<16) | (rand()& 0xffff) ; // reaplce this rand with something that has enough bits!!
  //  cout<<"size of int"<<sizeof(int) <<"\n";
// std::bitset<32>      x(k);
  
  //  std::cout << x << std::endl;

  /*  cout<<"bits of p\n";
    
    std::bitset<32>      pp(p->params[i]);
    std::cout << pp << std::endl;
        

    cout<<"bits of q\n";
    
    std::bitset<32>      qq(q->params[i]);
    std::cout << qq << std::endl;
        */

    replace->params[i] = (p->params[i]&k) | (q->params[i]& (~k)) ;
    
  // for(unsigned b)
 
   /* cout<<"new params" <<"\n";
        
    std::bitset<32>      pa(replace->params[i]);
    std::cout << pa<< std::endl;
    */
    
  }
    
}
void mutation(individual *p , double mutationrate){
  for(unsigned i =0; i< p->params.size() ;++i){
    for(unsigned b=0;b<sizeof(unsigned);++b){
      if( (rand() / double(RAND_MAX)) < mutationrate){
        p->params[i] = p->params[i] ^ (1<<b);
      }
    }
  /*  if(rand() / double(RAND_MAX) < mutationrate){
        int mask = (1<<  (rand() % (sizeof(int)*8 )) );
        p->params[i] = p->params[i] ^ mask;
      }
    */
   
  }
}
void reproduce(std::vector<individual*> &selected, std::vector<individual*> &replace, double crossrate,double mutationrate, std::uniform_int_distribution<int> &uni,std::mt19937 &rng){

  for(unsigned i =0;i<selected.size();++i ){
    crossover( selected[uni(rng)], selected[uni(rng)] ,replace[i] , crossrate );
    mutation(replace[i],mutationrate);
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
 
// [[Rcpp::export]]
Rcpp::List my_nsga2(Rcpp::Function func,Rcpp::NumericVector params,Rcpp::NumericVector objectives ,
                    Rcpp::NumericVector generations,Rcpp::NumericVector popsize,Rcpp::NumericVector cros,
                    Rcpp::NumericVector mutation, Rcpp::NumericVector lower, Rcpp::NumericVector upper){
    
    // shoudl check input
    if( params.size() != 1 ||popsize.size() != 1 ){
      cout<<"No value for params set\n";
    }
    unsigned int psize = popsize[0];
    unsigned int ngenerations =generations[0];
    unsigned int nparams = params[0];
    unsigned int nobj = objectives[0];
    double lowerbound = lower[0]; // fix the bound
    double upperbound = upper[0];
    double crossrate = cros[0];
    double mutationrate = mutation[0];
    
    std::list<std::list<individual*> > fronts ;
       
    std::vector<individual> population(psize*2) ;
    
 //   cout<<"start here\n";
    
    std::vector<individual*> parents(psize); // these to I need
    std::vector<individual*> replace(psize);    
    
    initPopulation(nparams, nobj ,population,parents,replace); // allocates memory
  //  cout<<"population inited\n";

    
     // init this one maybe in initpolupation
    //for( unsigned i=0;i<psize;++i){
  //    replace[i]= &population[i];
  //  }
  
    evaluate(parents, nparams,func, lowerbound, upperbound);
  //       cout<<"evaluate!\n";

    fastNonDominatedSort(population,0, psize,fronts);
 //   cout<<"fast non dom sort!\n";
    // select
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,psize-1); // guaranteed unbiased
 
    std::vector<individual*> selected(psize);

    for( unsigned i=0;i<psize;++i ){// should be vectorized somehow
      unsigned a = uni(rng);
      unsigned b = uni(rng);
      if(population[a].rank < population[b].rank){
        selected[i] = &population[a];
      }else{
        selected[i] = &population[b];
      }
    }
   // cout<<"tournament selection\n";

 //   for( unsigned i =psize;i<psize*2;++i){ // init this in init population too.
  //    replace[i-psize]=&population[i];
//    }
  //  cout<<"setup replace list\n";


        
    reproduce(selected,replace, crossrate,mutationrate,uni,rng);
    evaluate(replace, nparams,func, lowerbound, upperbound);
       
 //   cout<<"reproduce\n";
    


    unsigned P2 =  psize*2;
    for(unsigned i =0;i<ngenerations;++i){
      //population is already union of both...
       fastNonDominatedSort(population,0,P2,fronts);
       
  //     cout<<"frons estimated :"<<fronts.size()<<"\n";
       select_parents(fronts,nobj,parents ,replace);
       
    //cout<<"select parents\n";
  
        //tournement selection sould be vectorized somhow
       for( unsigned i=0;i<psize;++i ){
          unsigned a = uni(rng);
          unsigned b = uni(rng);
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
      //  cout<<"selected:end\n"<<"\n";
                    
       reproduce(selected,replace, crossrate,mutationrate,uni,rng);
       evaluate(replace, nparams,func, lowerbound, upperbound);
       
    }
    fastNonDominatedSort(population,0,P2,fronts);
    select_parents(fronts,nobj,parents ,replace); 
    // parents has now the results.
//    for(unsigned i =0;i<parents.size();++i){
//      cout<<"lokk rank:"<<parents[i]->rank<<"\n";
//    }
    
    Rcpp::NumericMatrix numerr(psize,nobj );
    Rcpp::NumericMatrix numparams(psize,nparams );
    
    std::vector<double> decparams(nparams);
    std::vector<int> frank(psize);
    
    
    //prepare results.
    for(unsigned k=0;k<psize;++k){
      decode( parents[k]->params , lowerbound, upperbound,decparams);
      frank[k] = parents[k]->rank;
   //   cout<<"rank "<<parents[k]->rank<<"\n";
      
      for(unsigned a=0;a<nparams;++a){
        numparams(k,a) = decparams[a];
      }    
      for(unsigned a=0;a<nobj;++a){
        numerr(k,a) = parents[k]->err[a];
      }    
      
    }
    
    //deletePolpulation(population); // free memory
    return Rcpp::List::create(Rcpp::Named("error") = numerr,
                              Rcpp::Named("parameters") = numparams,
                              Rcpp::Named("rank") = frank
                               
                              );
}


/// create multi treaded version
// [[Rcpp::export]]
Rcpp::List my_nsga2MT(Rcpp::Function func,Rcpp::NumericVector params,Rcpp::NumericVector objectives ,
                    Rcpp::NumericVector generations,Rcpp::NumericVector popsize,Rcpp::NumericVector cros,
                    Rcpp::NumericVector mutation, Rcpp::NumericVector lower, Rcpp::NumericVector upper,Rcpp::NumericVector Nthreads){
    
    // shoudl check input
    if( params.size() != 1 ||popsize.size() != 1 ){
      cout<<"No value for params set\n";
    }
    unsigned int psize = popsize[0];
    unsigned int ngenerations =generations[0];
    unsigned int nparams = params[0];
    unsigned int nobj = objectives[0];
    double lowerbound = lower[0]; // fix the bound
    double upperbound = upper[0];
    double crossrate = cros[0];
    double mutationrate = mutation[0];
    
    
    std::vector<individual> population(psize*2) ;
            
    /// set up threads stuff
    
     unsigned numThreadsX = Nthreads[0];
     if(numThreadsX == 0){
       numThreadsX =1;
     }
      if(numThreadsX > psize){
        cout<<"numThreads > psize setting numThreads = C\n";
        numThreadsX = psize;
      }
  
      unsigned int each = psize/numThreadsX;
                
      cout<<"threads:"<<numThreadsX<<" each has to evaluate:" <<each <<"\n";
      std::vector<pthread_t> threads(numThreadsX);
      std::vector<evalThread_data> tdata(numThreadsX);
      
    
      // set up data stucture
       for(unsigned int c=0;c<numThreadsX;++c){
       // tdata[c].population= &population;
        tdata[c].lower = lowerbound;
        tdata[c].uper = upperbound;
        tdata[c].nparams=nparams;
          
      //  rdata[c].decodeprams.resize(nparams);
         tdata[c].func = &func;
   
        
          tdata[c].from=c*each;
          tdata[c].to=c*each + each;
          if(c== numThreadsX-1){
             tdata[c].to=psize;
          }
       }
      pthread_attr_t attr;

      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    ///
    
    std::list<std::list<individual*> > fronts ;
    
 //   cout<<"start here\n";
    
    std::vector<individual*> parents(psize); // these to I need
    std::vector<individual*> replace(psize);    
    
    initPopulation(nparams, nobj ,population,parents,replace); // allocates memory
  //  cout<<"population inited\n";
    
     // init this one maybe in initpolupation
    //for( unsigned i=0;i<psize;++i){
  //    replace[i]= &population[i];
  //  }
  
 //   evaluate(parents, nparams,func, lowerbound, upperbound);
 
     cout<<"call evaluate!\n";
    
      evaluateMT(&parents,tdata,threads,&attr);
   
       cout<<"evaluate done\n";

    fastNonDominatedSort(population,0, psize,fronts);
 //   cout<<"fast non dom sort!\n";
    // select
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,psize-1); // guaranteed unbiased
 
    std::vector<individual*> selected(psize);

    for( unsigned i=0;i<psize;++i ){// should be vectorized somehow
      unsigned a = uni(rng);
      unsigned b = uni(rng);
      if(population[a].rank < population[b].rank){
        selected[i] = &population[a];
      }else{
        selected[i] = &population[b];
      }
    }
   // cout<<"tournament selection\n";

 //   for( unsigned i =psize;i<psize*2;++i){ // init this in init population too.
  //    replace[i-psize]=&population[i];
//    }
  //  cout<<"setup replace list\n";


        
    reproduce(selected,replace, crossrate,mutationrate,uni,rng);
    
    evaluate(replace, nparams,func, lowerbound, upperbound);
       
 //   cout<<"reproduce\n";
    


    unsigned P2 =  psize*2;
    for(unsigned i =0;i<ngenerations;++i){
      //population is already union of both...
       fastNonDominatedSort(population,0,P2,fronts);
       
  //     cout<<"frons estimated :"<<fronts.size()<<"\n";
       select_parents(fronts,nobj,parents ,replace);
       
    //cout<<"select parents\n";
  
        //tournement selection sould be vectorized somhow
       for( unsigned i=0;i<psize;++i ){
          unsigned a = uni(rng);
          unsigned b = uni(rng);
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
      //  cout<<"selected:end\n"<<"\n";
                    
       reproduce(selected,replace, crossrate,mutationrate,uni,rng);
       evaluate(replace, nparams,func, lowerbound, upperbound);
       
    }
    fastNonDominatedSort(population,0,P2,fronts);
    select_parents(fronts,nobj,parents ,replace); 
    // parents has now the results.

    Rcpp::NumericMatrix numerr(psize,nobj );
    Rcpp::NumericMatrix numparams(psize,nparams );
    
    std::vector<double> decparams(nparams);
    std::vector<int> frank(psize);
    
    
    //prepare results.
    for(unsigned k=0;k<psize;++k){
      decode( parents[k]->params , lowerbound, upperbound,decparams);
      frank[k] = parents[k]->rank;
      //   cout<<"rank "<<parents[k]->rank<<"\n";
      
      for(unsigned a=0;a<nparams;++a){
        numparams(k,a) = decparams[a];
      }    
      for(unsigned a=0;a<nobj;++a){
        numerr(k,a) = parents[k]->err[a];
      }    
      
    }
    
    //deletePolpulation(population); // free memory
    return Rcpp::List::create(Rcpp::Named("error") = numerr,
                              Rcpp::Named("parameters") = numparams,
                              Rcpp::Named("rank") = frank
                               
                              );
}


////////////////    


////try vectorized!!



// [[Rcpp::export]]
Rcpp::List my_nsga2Vec(Rcpp::Function func,Rcpp::NumericVector params,Rcpp::NumericVector objectives ,
                    Rcpp::NumericVector generations,Rcpp::NumericVector popsize,Rcpp::NumericVector cros,
                    Rcpp::NumericVector mutation, Rcpp::NumericVector lower, Rcpp::NumericVector upper,Rcpp::Function repfunc){
    
    // shoudl check input
    if( params.size() != 1 ||popsize.size() != 1 ){
      cout<<"No value for params set\n";
    }
    srand (time(NULL)); // init rnd
    
    unsigned int psize = popsize[0];
    unsigned int ngenerations =generations[0];
    unsigned int nparams = params[0];
    unsigned int nobj = objectives[0];
    double lowerbound = lower[0]; // fix the bound
    double upperbound = upper[0];
    double crossrate = cros[0];
    double mutationrate = mutation[0];
    
    std::list<std::list<individual*> > fronts ;
       
    std::vector<individual> population(psize*2) ;
    
 //   cout<<"start here\n";
    
    std::vector<individual*> parents(psize); // these to I need
    std::vector<individual*> replace(psize);    
    
    initPopulation(nparams, nobj ,population,parents,replace); // allocates memory
  //  cout<<"population inited\n";

    evaluateVec(parents, nparams,func, lowerbound, upperbound,nobj);
  //       cout<<"evaluate!\n";

    fastNonDominatedSort(population,0, psize,fronts);
 //   cout<<"fast non dom sort!\n";
    // select
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,psize-1); // guaranteed unbiased
 
    std::vector<individual*> selected(psize);

    for( unsigned i=0;i<psize;++i ){// should be vectorized somehow
      unsigned a = uni(rng);
      unsigned b = uni(rng);
      if(population[a].rank < population[b].rank){
        selected[i] = &population[a];
      }else{
        selected[i] = &population[b];
      }
    }
   // cout<<"tournament selection\n";

 //   for( unsigned i =psize;i<psize*2;++i){ // init this in init population too.
  //    replace[i-psize]=&population[i];
//    }
  //  cout<<"setup replace list\n";


        
    reproduce(selected,replace, crossrate,mutationrate,uni,rng);
    evaluateVec(replace, nparams,func, lowerbound, upperbound,nobj);
       
 //   cout<<"reproduce\n";
    


    unsigned P2 =  psize*2;
    for(unsigned i =0;i<ngenerations;++i){
      //population is already union of both...
       fastNonDominatedSort(population,0,P2,fronts);
       
  //     cout<<"frons estimated :"<<fronts.size()<<"\n";
       select_parents(fronts,nobj,parents ,replace);
       
       reportCurrentPopErr(parents,nobj,  repfunc);
       
    //cout<<"select parents\n";
  
        //tournement selection sould be vectorized somhow
       for( unsigned i=0;i<psize;++i ){
          unsigned a = uni(rng);
          unsigned b = uni(rng);
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
      //  cout<<"selected:end\n"<<"\n";
                    
       reproduce(selected,replace, crossrate,mutationrate,uni,rng);
       evaluateVec(replace, nparams,func, lowerbound, upperbound,nobj);
       
    }
    fastNonDominatedSort(population,0,P2,fronts);
    select_parents(fronts,nobj,parents ,replace); 
    
    reportCurrentPopErr(parents,nobj,  repfunc);
    
    // parents has now the results.
//    for(unsigned i =0;i<parents.size();++i){
//      cout<<"lokk rank:"<<parents[i]->rank<<"\n";
//    }
    
    Rcpp::NumericMatrix numerr(psize,nobj );
    Rcpp::NumericMatrix numparams(psize,nparams );
    
    std::vector<double> decparams(nparams);
    std::vector<int> frank(psize);
    
    
    //prepare results.
    for(unsigned k=0;k<psize;++k){
      decode( parents[k]->params , lowerbound, upperbound,decparams);
      frank[k] = parents[k]->rank;
   //   cout<<"rank "<<parents[k]->rank<<"\n";
      
      for(unsigned a=0;a<nparams;++a){
        numparams(k,a) = decparams[a];
      }    
      for(unsigned a=0;a<nobj;++a){
        numerr(k,a) = parents[k]->err[a];
      }
      
    }
    
    //deletePolpulation(population); // free memory
    return Rcpp::List::create(Rcpp::Named("error") = numerr,
                              Rcpp::Named("parameters") = numparams,
                              Rcpp::Named("rank") = frank
                               
                              );
}
