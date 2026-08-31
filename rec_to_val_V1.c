#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>

#include <omp.h>
#include <string.h>
#include <flint/mpn_extras.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/arf.h>
#include <flint/profiler.h>
#include <flint/thread_support.h>

typedef struct
{
	fmpz_poly_mat_t * Rec; slong r_lcoeff; fmpz_mat_t * Num_inhom; fmpz* Den_inhom; fmpz_mat_t * Num_inits; fmpz* Den_inits; 
		bool have_inhom; slong next_prime; slong CHUNK; slong PACKED_CHUNK; slong BITS;  slong r_ord; slong rec_degree; slong n_val; slong start_init; 
		slong start_comp; slong length_init; slong end; slong k_write; nn_srcptr primes; nmod_t * nmod_primes; const char* filename_template_temp; slong num_worker;
}
threaded_args;


void compute_rec_values_threaded(slong i,void * args);

__attribute__ ((hot))
void compute_rec_values(fmpz_poly_mat_t  Rec,slong r_lcoeff,fmpz_mat_t * Num_inhom,fmpz* Den_inhom,fmpz_mat_t * Num_inits,fmpz* Den_inits,
		bool have_inhom,slong next_prime,slong CHUNK,slong PACKED_CHUNK,slong BITS, slong r_ord,slong rec_degree,slong n_val,slong start_init,
		slong start_comp,slong length_init,slong end,slong k_write,nn_srcptr primes,nmod_t * nmod_primes,const char* filename_template_temp,slong num_worker);


ulong poly_evaluate_nmod_powprecomp(const nmod_poly_t poly,nn_srcptr c_pow,dot_params_t dotpar);
slong readFilePoly(fmpz_poly_mat_t Rec,slong* coeff,const char * filename_rec);
slong readFileRationals(fmpz_mat_t * Val_Num,fmpz** Denominators,slong n_val, slong end,const char * filename_initial);
char * get_temp_file_name(char * output,const char * filename_template_temp,slong i);


int main(int argc, char* argv[]){
omp_get_num_procs();
const slong CHUNK=64;
//const ulong firstprime=1512762481;//1514687359;//< 2^(30.5)
//const ulong BITS=31;
const ulong firstprime=6085156183315683391; //<2^(62.5)
const ulong BITS=63;
const ulong PACKED_CHUNK=(CHUNK * BITS - 1) / FLINT_BITS + 1;
//const ulong firstprime=18446744073708659869U; //<2^(64)

//const slong MAX_NUM_PRIMES=32768;//If more primes are needed this might point to an error


//flint_printf("num_threads: %wd\n",flint_get_num_threads());

//slong start;//index of the first initial value
slong k; //index of the next value which will be computed
slong i,j,p,l;
slong next_prime=0;
bool increase_primes=true;


char * filename_out;
char * filename_rec;
char * filename_deq;
char * filename_template_temp;
char filename_temp[300];
char * filename_initial;
char * filename_inhom;
FILE * file_out;
if(argc != 10 && argc != 11 && argc != 12){
	printf("Usage: %s <filename_out> <filename_rec> <filename_deq> <filename_template_temp> <n_val> <start_comp> <end> <final_precision> <num_threads>\n", argv[0]);
    printf("   Or: %s <filename_out> <filename_rec> <filename_deq> <filename_template_temp> <n_val> <start_comp> <end> <final_precision> <num_threads> <filename_initial>\n", argv[0]);
    printf("   Or: %s <filename_out> <filename_rec> <filename_deq> <filename_template_temp> <n_val> <start_comp> <end> <final_precision> <num_threads> <filename_initial> <filename_inhom>\n", argv[0]);
    return 12;
}
filename_out = argv[1];
filename_rec = argv[2];
filename_deq = argv[3];
filename_template_temp = argv[4];
const slong n_val = atol(argv[5]);  //number of sets of initial values
const slong start_comp = atol(argv[6]);
const slong end = atol(argv[7]);   //last value which should be computed
const slong final_precision = atol(argv[8]); 
const slong num_threads = atol(argv[9]); 
slong n_primes=num_threads*CHUNK; //initial number of homomorphic images. n_primes%CHUNK==0
FILE ** file_tmp=flint_malloc(n_primes*sizeof(FILE *));
if (argc >= 11){
    filename_initial = argv[10];
    flint_printf("filename_initial: %s\n", filename_initial);
}
if (argc >= 12){
    filename_inhom= argv[11]; 
    flint_printf("filename_inhom: %s\n", filename_inhom);
}
const bool have_init_vals=(argc >=11); 
const bool have_inhom=(argc >= 12); 

flint_printf("Computed with FLINT-%s\n", flint_version);
# if defined(__AVX__)
    printf("Using AVX2\n");
# endif
flint_printf("filename_out: %s\n", filename_out);
flint_printf("filename_rec: %s\n", filename_rec);
flint_printf("filename_deq: %s\n", filename_deq);
flint_printf("filename_template_temp: %s\n", filename_template_temp);
flint_printf("n_val: %wd\n", n_val);
flint_printf("start_comp: %wd\n", start_comp);
flint_printf("end: %wd\n", end);
flint_printf("final_precision(base 2): %wd\n", final_precision);
flint_printf("num_threads %wd\n", num_threads);

//Avoid overwriting filename_out
if (access(filename_out, F_OK) == 0) {
    flint_printf("output file exists already: %s\n",filename_out);
    return 13;
}

//Setting filenames of temporary files
get_temp_file_name(filename_temp,filename_template_temp,0);
printf("first temp-file: %s\n",filename_temp);


fmpz_poly_mat_t Rec;
//fmpz_poly_mat_t Deq;


fmpz_t output,input,den,den_add,num_bound;
//current denominator for reconstruction
fmpz_t modulus;//product of all primes
fmpz_mat_t Num_inits,CurValsZ;
fmpz * Den_inits;
fmpq_mat_t CurValsQ;
fmpz_comb_t comb;
fmpz_comb_temp_t comb_temp;


arf_t out_numeric;
fmpz_init(output);fmpz_init(input);fmpz_init(den);fmpz_init(num_bound);fmpz_init(den_add);
arf_init(out_numeric);
fmpz_init_set_si(modulus,WORD(1));

//nmod_mat_t Values_k_mod[MAX_NUM_PRIMES];
nmod_mat_t * Values_k_mod=flint_malloc(n_primes*sizeof(nmod_mat_t));
nn_ptr val_chunk_k=_nmod_vec_init(CHUNK);
nn_ptr val_chunk_packed=_nmod_vec_init(PACKED_CHUNK);
nn_ptr den_mod=_nmod_vec_init(n_primes);


nn_ptr primes=_nmod_vec_init(n_primes);
nmod_t* nmod_primes= flint_malloc(n_primes*sizeof(nmod_t));
primes[0] = firstprime;
nmod_init(&nmod_primes[0],firstprime);
for (p = 1; p < n_primes; p++){//maybe use n_randprime?
    primes[p] = n_nextprime(primes[p-1], 0);
    nmod_init(&nmod_primes[p],primes[p]);
}

TIMEIT_ONCE_START
slong r_lcoeff;
//slong d_tcoeff;
const slong r_ord=readFilePoly(Rec,&r_lcoeff,filename_rec);
if(r_ord<0)
    return r_ord;
fmpz_poly_neg(fmpz_poly_mat_entry(Rec, 0, r_ord),fmpz_poly_mat_entry(Rec, 0, r_ord));


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                          set initial values
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slong start_init,length_init;
if(have_init_vals){
    start_init=readFileRationals(&Num_inits,&Den_inits,n_val,start_comp,filename_initial);
    length_init=start_comp-start_init;
    //SHOW_MEMORY_USAGE;
    print_memory_usage();
} else {
	//fmpz_poly_mat_t Deq;
	slong d_tcoeff=start_comp-n_val;
	start_init=0;	
	slong d_ord=start_comp;
    //slong d_ord=readFilePoly(Deq,&d_tcoeff,filename_deq);
    //fmpz_poly_neg(fmpz_poly_mat_entry(Deq, 0, d_ord),fmpz_poly_mat_entry(Deq, 0, d_ord));
    length_init=start_comp-start_init;
    //if(d_ord<0)
        //return d_ord; 
    fmpz_mat_init(Num_inits, n_val, length_init);
    fmpz_mat_zero(Num_inits);
    Den_inits=_fmpz_vec_init(length_init);
    for(slong i=0;i<length_init;i++)
		fmpz_set_si(Den_inits+i,WORD(1));
    for(i=d_tcoeff;i<n_val+d_tcoeff;i++){
        for(j=0;j<d_ord;j++){
            if(j==i)
                fmpz_set_si(fmpz_mat_entry(Num_inits,i-d_tcoeff,j),WORD(1));
        }
    }
    flint_printf("d_ord: %wd \n",d_ord);
}
fmpz_set(den,Den_inits+length_init-1); 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                          set inhom
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmpz_mat_t Num_inhom;
fmpz * Den_inhom;

if(have_inhom){
	readFileRationals(&Num_inhom,&Den_inhom,n_val,end+1,filename_inhom);
	//SHOW_MEMORY_USAGE
	print_memory_usage();
}





//const slong deq_degree=fmpz_poly_mat_max_length(Deq);//Actually degree+1
const slong rec_degree=fmpz_poly_mat_max_length(Rec);//Actually degree+1
slong k_write=start_comp;    //next value for which we have to do reconstruction

/*
nmod_mat_t Values_mod[CHUNK];
nn_ptr Rec_mod_k=_nmod_vec_init(r_ord+1);
nn_ptr Next_val_mod=_nmod_vec_init(n_val);
nmod_poly_mat_t Rec_mod;
nmod_mat_t Values_mod_win,Rec_mod_mat;
nn_ptr k_pow=_nmod_vec_init(rec_degree);
nn_ptr num_inhom_k_mod=_nmod_vec_init(n_val);
//nmod_mat_t K_fpow;
//nmod_mat_init(K_fpow,start+r_ord+1,d_ord+1,primes[0]);
for(l=0;l<CHUNK;l++)
	nmod_mat_init(Values_mod[l],n_val,end-start_init+1,primes[0]);
SHOW_MEMORY_USAGE
*/

flint_printf("r_ord: %wd \n",r_ord);
flint_printf("r_lcoeff: %wd \n",r_lcoeff);
//flint_printf("d_tcoeff: %wd \n",d_tcoeff);
//flint_printf("deq_degree: %wd \n",deq_degree);
flint_printf("start_init: %wd \n",start_init);
flint_printf("start_comp: %wd \n",start_comp);
flint_printf("end: %wd \n",end);
flint_printf("n_val: %wd \n",n_val);
flint_printf("deg: %wd \n",rec_degree);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                          main loop 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(increase_primes){
    TIMEIT_ONCE_START  
    threaded_args  args[1];
	args->Rec = &Rec;
	args->r_lcoeff = r_lcoeff;
	args->Num_inhom = &Num_inhom;
	args->Den_inhom = Den_inhom;
	args->Num_inits = &Num_inits;
	args->Den_inits = Den_inits;
	args->have_inhom = have_inhom;
	args->next_prime = next_prime;
	args->CHUNK = CHUNK;
	args->PACKED_CHUNK = PACKED_CHUNK;
	args->BITS = BITS;
	args->r_ord = r_ord;
	args->rec_degree = rec_degree;
	args->n_val = n_val;
	args->start_init = start_init;
	args->start_comp = start_comp;
	args->length_init = length_init;
	args->end = end;
	args->k_write = k_write;
	args->primes = primes;
	args->nmod_primes = nmod_primes;
	args->filename_template_temp = filename_template_temp;
	
	flint_set_num_threads(num_threads);
	flint_parallel_do(compute_rec_values_threaded,args,num_threads,num_threads,FLINT_PARALLEL_UNIFORM);
	flint_set_num_threads(1);


//	compute_rec_values(Rec, r_lcoeff,& Num_inhom,Den_inhom,&Num_inits,Den_inits,
	//	have_inhom,next_prime,CHUNK,PACKED_CHUNK,BITS,r_ord,rec_degree,n_val, start_init,
	//	 start_comp, length_init, end, k_write, primes,nmod_primes,filename_template_temp,0);
    
/*    for(p=next_prime;p<n_primes;p++){//this loop may be parallized 
        //Set modulus to current prime
        nmod_t mod=nmod_primes[p];
        nmod_mat_set_mod(Values_mod[p%CHUNK],primes[p]);
        nmod_mat_init(Rec_mod_mat,r_ord+1,rec_degree,primes[p]);
        //nmod_mat_set_mod(Rec_mod_mat,primes[p]);
        
        nmod_poly_mat_init(Rec_mod, 1, r_ord+1,primes[p]); 
        const dot_params_t rec_eval_dotpar=_nmod_vec_dot_params(rec_degree,mod);
        
        TIMEIT_ONCE_START
        
        
        //reduce Rec modulo mod
        for(i=0;i<=r_ord;i++)
            fmpz_poly_get_nmod_poly(nmod_poly_mat_entry(Rec_mod,0,i),fmpz_poly_mat_entry(Rec,0,i));

		for(i=0;i<=r_ord;i++)
			for(j=0;j<rec_degree;j++)
				nmod_mat_entry(Rec_mod_mat,i,j)=nmod_poly_get_coeff_ui(nmod_poly_mat_entry(Rec_mod,0,i),j);
				
        //Set Values_mod to Values modulo mod
        nmod_mat_window_init(Values_mod_win,Values_mod[p%CHUNK],0,0,n_val,length_init);
        fmpz_mat_get_nmod_mat(Values_mod_win,Num_inits);
        for(slong k=0;k<length_init;k++){        
			ulong den_inv_mod=nmod_inv(fmpz_get_nmod(Den_inits+k,mod),mod);
			for(slong l=0;l<n_val;l++){
				nmod_mat_entry(Values_mod_win,l,k)=nmod_mul(nmod_mat_entry(Values_mod_win,l,k),den_inv_mod,mod);
			}
        }
        nmod_mat_window_clear(Values_mod_win);
        
        ////multipoint eval
        //nn_ptr * tree=_nmod_poly_tree_alloc(end-start_init+1);
        //nn_ptr range=_nmod_vec_init(end-start_init+1);
        //nmod_mat_t Evals;
        //nmod_mat_init(Evals,r_ord+1,end-start_init+1,mod.n);
        //for(slong k=start_init;k<=end;k++)
			//range[k-start_init]=k;//reduced!!!!
        //_nmod_poly_tree_build(tree,range,end-start_init+1,mod);
        //for(slong i=0;i<=r_ord;i++)
			//_nmod_poly_evaluate_nmod_vec_fast_precomp(&nmod_mat_entry(Evals,i,0),nmod_poly_mat_entry(Rec_mod,0,i)->coeffs,rec_degree,tree,end-start_init+1,mod);
        //_nmod_poly_tree_free(tree,end-start_init+1);
        //nmod_mat_clear(Evals);
        //_nmod_vec_clear(range);
        //------------------------------------
        //             Rec loop
        //------------------------------------
        TIMEIT_ONCE_START
        k_pow[0]=1;
        for(k = start_comp; k <= end; k++){
            for(i=1;i<rec_degree;i++)
                k_pow[i]=nmod_mul(k_pow[i-1],k-r_lcoeff,mod);
            slong n_previous=FLINT_MIN(r_ord,k-start_init);
            for(i=r_ord-n_previous;i<=r_ord;i++){ //Rec_mod_k=Rec_mod/.n->k
				Rec_mod_k[i]=_nmod_vec_dot(&nmod_mat_entry(Rec_mod_mat,i,0),k_pow,rec_degree,mod,rec_eval_dotpar);
                //Rec_mod_k[i]=poly_evaluate_nmod_powprecomp(nmod_poly_mat_entry(Rec_mod,0,i),k_pow,rec_eval_dotpar);
                //Rec_mod_k[i]=nmod_poly_evaluate_nmod(nmod_poly_mat_entry(Rec_mod,0,i),k); 
            }
            nmod_mat_window_init(Values_mod_win,Values_mod[p%CHUNK],0,k-start_init-n_previous,n_val,k-start_init);
            nmod_mat_mul_nmod_vec(Next_val_mod,Values_mod_win,Rec_mod_k+r_ord-n_previous,n_previous);
            if(have_inhom){
				 _fmpz_vec_get_nmod_vec(num_inhom_k_mod,fmpz_mat_entry(Num_inhom,k-start_comp,0),n_val,mod); 
				 ulong inverse= nmod_inv(fmpz_get_nmod(Den_inhom+k-start_comp,mod),mod);
				 _nmod_vec_scalar_mul_nmod(num_inhom_k_mod,num_inhom_k_mod,n_val,inverse,mod);
				 _nmod_vec_add(Next_val_mod,Next_val_mod,num_inhom_k_mod,n_val,mod);
			}
            if(have_inhom){
                for(slong j=0;j<n_val;j++)
				    num_inhom_k_mod[j]=fmpz_get_nmod(fmpz_mat_entry(Num_inhom,j,k-start_comp),mod); 
				 ulong inverse= nmod_inv(fmpz_get_nmod(Den_inhom+k-start_comp,mod),mod);
				 _nmod_vec_scalar_mul_nmod(num_inhom_k_mod,num_inhom_k_mod,n_val,inverse,mod);
				 _nmod_vec_add(Next_val_mod,Next_val_mod,num_inhom_k_mod,n_val,mod);
			}
            
            
            for(j=0;j<n_val;j++){
                nmod_mat_entry(Values_mod[p%CHUNK],j,k-start_init)=nmod_div(Next_val_mod[j],Rec_mod_k[r_ord],mod);
            } 
            nmod_mat_window_clear(Values_mod_win);
        } 
        nmod_mat_clear(Rec_mod_mat);
        nmod_poly_mat_clear(Rec_mod); 
        printf("  loop rec: ");
        TIMEIT_ONCE_STOP
        
        
        //write computed residues to file. Values_mod will be overwritten later to save memory
        if((p+1)%CHUNK==0){
            TIMEIT_ONCE_START               
            FILE * file_tmp = fopen(get_temp_file_name(filename_temp,filename_template_temp,p/CHUNK), "w");
            for(k = k_write; k <= end; k++){
                 //flint_printf("write k = %wd\n",k);
                for(j=0;j<n_val;j++){ 
                    for(l=0;l<CHUNK;l++){
                        val_chunk_k[l]=nmod_mat_get_entry(Values_mod[l],j,k-start_init);
                    }                    
                    _nmod_poly_bit_pack(val_chunk_packed, val_chunk_k, CHUNK, BITS);
                    fmpz_set_ui_array(output,val_chunk_packed,PACKED_CHUNK);
                    //fmpz_set_ui_array(output,val_chunk_k,CHUNK);
                    if(fmpz_out_raw(file_tmp,output)==0){
                        flint_printf("unable to write to file\n");
                        return 15;
                    }
                }       
            }
            
            
            fclose(file_tmp);
            printf("  write to file: ");
            TIMEIT_ONCE_STOP
        }
        printf("computed residue %lu: ",p);
        TIMEIT_ONCE_STOP
        if(end>3260-start_init && n_val>1)
		    flint_printf("Value: %wu \n",nmod_mat_get_entry(Values_mod[p%CHUNK],1,3260-start_init));
		
    }*/

    printf("Computing residues combined: ");
    TIMEIT_ONCE_STOP;
    //SHOW_MEMORY_USAGE;
    print_memory_usage();


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//         do reconstruction
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    TIMEIT_ONCE_START
    for (p = FLINT_MAX(next_prime-1,0); p < n_primes-1; p++)
        fmpz_mul_ui(modulus,modulus,primes[p]);
    for(p=next_prime;p<n_primes;p++)
        nmod_mat_init(Values_k_mod[p],1,n_val,primes[p]);
    nmod_t mod_max;
    nmod_init(&mod_max,UWORD_MAX);
    fmpz_t mantissa,exponent;
	fmpz_init(mantissa);fmpz_init(exponent);
    /* Data needed by multi CRT functions */
    fmpz_comb_init(comb, primes, n_primes-1);
    fmpz_comb_temp_init(comb_temp, comb);
	fmpz_mat_init(CurValsZ,1,n_val);
	fmpq_mat_init(CurValsQ,1,n_val);
    for(l=next_prime/CHUNK;l<n_primes/CHUNK;l++)
        file_tmp[l]=fopen(get_temp_file_name(filename_temp, filename_template_temp,l) ,"r");
    file_out=fopen(filename_out,"a");
    //nmod_t checkmod;

    fmpz_t den_bound;fmpz_init(den_bound);
	
    fmpz_t den_bound_safety;fmpz_init(den_bound_safety);

    //nmod_init(&checkmod,primes[n_primes-1]);
    for(p=next_prime;p<n_primes;p++)
        den_mod[p]=fmpz_get_nmod(den,nmod_primes[p]);

    for(k = k_write;k <= end; k++){
        k_write=k;
        TIMEIT_ONCE_START
        //read homomorphic images from tmp-files
        l=(increase_primes?(next_prime/CHUNK):0);//we continue where we left before we computed more primes      
        for(;l<n_primes/CHUNK;l++){
            for(j=0;j<n_val;j++){
                if(fmpz_inp_raw(input,file_tmp[l])==0){
                    flint_printf("unable to read from file: %s\n",get_temp_file_name(filename_temp, filename_template_temp,l));
                    return 1;
                }
                //fmpz_get_ui_array(val_chunk_k,CHUNK,input);
                fmpz_get_ui_array(val_chunk_packed,PACKED_CHUNK,input);
                _nmod_poly_bit_unpack(val_chunk_k,CHUNK,val_chunk_packed,BITS,mod_max);
                for(i=0;i<CHUNK;i++){
                    nmod_mat_entry(Values_k_mod[l*CHUNK+i],0 ,j)=nmod_mul(val_chunk_k[i],den_mod[l*CHUNK+i],nmod_primes[l*CHUNK+i]);
                }
            }
        }
        //bound for the denominator
        if(k<length_init+200){//at the beginning (especially k=d_ord) the denominater might be unpredictable. Set large bound on den in this case
            fmpz_one_2exp(den_bound,fmpz_sizeinbase(modulus,2)/2-256);
            fmpz_set(num_bound,den_bound);        
        } else {
            fmpz_one_2exp(den_bound,512);//2^512, by which factor the denominator increases at most in each step (can be increasesed significantly if it is to low)
            fmpz_cdiv_q_2exp(num_bound,modulus,1024);// modulus/2^1024, this gives 50 binary digits additional safety for reconstruction
        }      

        //reconstruct
        increase_primes=false;       
        fmpz_mat_multi_CRT_ui_precomp(CurValsZ, Values_k_mod,n_primes-1,comb,comb_temp,0);
        for(j=0;j<n_val;j++){
            
            fmpq * CurValsQ_j=fmpq_mat_entry(CurValsQ,0,j);
            int suc;
			suc=fmpq_reconstruct_fmpz_2(CurValsQ_j,fmpz_mat_entry(CurValsZ,0,j),modulus,num_bound,den_bound);
            //using last prime for checking correctness
            ulong test;
            if(suc){
				test=nmod_div(fmpz_get_nmod(fmpq_numref(CurValsQ_j),nmod_primes[n_primes-1]),
                           fmpz_get_nmod(fmpq_denref(CurValsQ_j),nmod_primes[n_primes-1]),nmod_primes[n_primes-1]);
            }
                      
            if(!suc || (test!=nmod_mat_entry(Values_k_mod[n_primes-1],0 ,j))||n_primes<20){
                increase_primes=true;
                flint_printf("failed for k=%wd increase number of primes\n",k);
                break;
            }
        }
        if(!increase_primes){
            //update den
            fmpq_mat_get_fmpz_mat_matwise(CurValsZ,den_add,CurValsQ);
            fmpz_mul(den,den,den_add);
            for(p=0;p<n_primes;p++)
                den_mod[p]=nmod_mul(den_mod[p],fmpz_get_nmod(den_add,nmod_primes[p]),nmod_primes[p]);
            
           
			if(final_precision>0){
				for(j=0;j<n_val;j++){                
					arf_fmpz_div_fmpz(out_numeric, fmpz_mat_entry(CurValsZ,0,j),den,final_precision+1, ARF_RND_NEAR);
					arf_get_fmpz_2exp(mantissa,exponent,out_numeric);
					fmpz_fprint(file_out,mantissa);
					flint_fprintf(file_out," ");   
					fmpz_fprint(file_out,exponent);
				    //arf_fprint(file_out,out_numeric);
					flint_fprintf(file_out," ");                
				}
				flint_fprintf(file_out,"\n"); 
            } else {// if  final_precision<=0 then write exact values
			    //write to file "filename_out". First the n_val numerators, then their common denominator
			    for(j=0;j<n_val;j++){
				   fmpz_fprint(file_out,fmpz_mat_entry(CurValsZ,0,j));
				   flint_fprintf(file_out," ");
				}
				fmpz_fprint(file_out,den);
				flint_fprintf(file_out,"\n");
			}
            
            flint_printf("computed residue %wd: ",k);
        } else {
            flint_printf("Reconstruction for k = %wd failed, computing more images.\n",k);
            next_prime=n_primes;
            n_primes+=num_threads*CHUNK;
            //make vectors longer
            nn_ptr temp=_nmod_vec_init(n_primes);
             _nmod_vec_set(temp,primes,next_prime);
            _nmod_vec_clear(primes);
            primes=temp;
            
            temp=_nmod_vec_init(n_primes);
             _nmod_vec_set(temp,den_mod,next_prime);
            _nmod_vec_clear(den_mod);
            den_mod=temp;
            
            nmod_mat_t * Values_k_mod_temp=flint_malloc(n_primes*sizeof(nmod_mat_t));
            for(slong i=0;i<next_prime;i++)
				nmod_mat_swap(Values_k_mod_temp[i],Values_k_mod[i]);
			flint_free(Values_k_mod);
            Values_k_mod=Values_k_mod_temp;   
            
            nmod_t* temp_mod= flint_malloc(n_primes*sizeof(nmod_t));
            for(slong i=0;i<next_prime;i++)
				temp_mod[i]=nmod_primes[i];
			flint_free(nmod_primes);
            nmod_primes=temp_mod;          
           for (slong p = next_prime; p < n_primes; p++){
				primes[p] = n_nextprime(primes[p-1], 0);
				nmod_init(&nmod_primes[p],primes[p]);
			}
			FILE ** temp_file_tmp=flint_malloc(n_primes/CHUNK*sizeof(FILE *));
			for(slong i=0;i<next_prime/CHUNK;i++)
				temp_file_tmp[i]=file_tmp[i];
			flint_free(file_tmp);
			file_tmp=temp_file_tmp;
            k=end+1;
        }
        TIMEIT_ONCE_STOP;
    }
    fmpz_clear(den_bound);fmpz_clear(mantissa);fmpz_clear(exponent);
    fmpz_clear(den_bound_safety);
    fmpz_mat_clear(CurValsZ);
    fmpq_mat_clear(CurValsQ);
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(comb_temp);
    fclose(file_out);
    TIMEIT_ONCE_STOP;

}


//delete temporary files
for(l=0;l<n_primes/CHUNK;l++){
    fclose(file_tmp[l]);
    remove(get_temp_file_name(filename_temp,filename_template_temp,l));
}



printf("Everything combined: ");
TIMEIT_ONCE_STOP;
//SHOW_MEMORY_USAGE;
print_memory_usage();

fmpz_poly_mat_clear(Rec);
fmpz_mat_clear(Num_inits);

fmpz_clear(den);
return 0;
}


////////////////////////////////////////////////////////////////////////////////////
slong readFileRationals(fmpz_mat_t * Val_Num,fmpz** Denominators,slong n_val, slong end,const char * filename_initial){
    printf("set initial values: \n");
    slong start;
    fmpq_mat_t CurVals;
    fmpz_mat_t Val_Num_i;
    FILE * file_initial = fopen(filename_initial, "r");
    if (file_initial == NULL){
        flint_printf("unable to open file of initial values %s\n",filename_initial);
        return -16;
    }
    flint_fscanf(file_initial, "%wd", &start);
    slong length=end - start;
    fmpq_mat_init(CurVals,n_val,1);
    fmpz_mat_init(*Val_Num,n_val,length);
    *Denominators=_fmpz_vec_init(length);
    //flint_printf("start_init %wd\n",start);
    //flint_printf("start_comp %wd\n",end);
    //flint_printf("length %wd\n",length);
    //flint_printf("n_val %wd\n",n_val);
    for(slong i=0;i<length;i++){
        for(slong j=0;j<n_val;j++){
            if(feof(file_initial)){  
                flint_printf("unable to read full file of initial values\n");
                return -17;
            }
            fmpz_fread(file_initial,fmpq_numref(fmpq_mat_entry(CurVals,j,0)));
            if(feof (file_initial)){  
                flint_printf("unable to read full file of initial values\n");
                return -17;
            }
            fmpz_fread(file_initial,fmpq_denref(fmpq_mat_entry(CurVals,j,0))); 
            fmpz_mat_window_init(Val_Num_i,*Val_Num,0,i,n_val,i+1);
            fmpq_mat_get_fmpz_mat_matwise(Val_Num_i,*Denominators+i,CurVals);
            fmpz_mat_window_clear(Val_Num_i);
        }
    }
    fclose(file_initial);
    
    fmpq_mat_clear(CurVals);
    return start;
}

slong readFilePoly(fmpz_poly_mat_t Rec,slong* coeff,const char * filename_rec){
    printf("set entries of recursion: \n");    
    //set entries of recursion
    slong i,ord;
    FILE * file_rec = fopen(filename_rec, "r");
    if (file_rec == NULL){
        flint_printf("unable to open file\n");
        return -1;
    }
    flint_fscanf(file_rec, "%wd", &ord);  
    flint_fscanf(file_rec, "%wd", coeff);  
    fmpz_poly_mat_init(Rec, 1, ord+1);
    for(i=0;i<=ord;i++){
        if(feof (file_rec)){  
            flint_printf("unable to read full file\n");
            return -1;
        }
        if(fmpz_poly_fread(file_rec,fmpz_poly_mat_entry(Rec, 0, i))==0){
            printf("unable to read full file\n");
            return -1;  
        }
    }
    fclose(file_rec);
    return ord;
}

ulong poly_evaluate_nmod_powprecomp(const nmod_poly_t poly,nn_srcptr c_pow,dot_params_t dotpar)
{
    if (poly->length == 0)
        return 0;
    return _nmod_vec_dot(poly->coeffs,c_pow,poly->length,poly->mod,dotpar);
}

char * get_temp_file_name(char * output,const char * filename_template_temp,slong i){
    strcpy(output,filename_template_temp); 
    strcat(output, "_");
    char temp[21];
    sprintf(temp, "%ld", i);
    strcat(output,temp);
    strcat(output,".tmp"); 
    return output;
}


void compute_rec_values_threaded(slong num_thread,void * args_ptr){
	threaded_args args = *((threaded_args *) args_ptr);
	compute_rec_values(*(args.Rec), args.r_lcoeff,args.Num_inhom,args.Den_inhom,args.Num_inits,args.Den_inits,
		args.have_inhom,args.next_prime,args.CHUNK,args.PACKED_CHUNK,args.BITS,args.r_ord,args.rec_degree,args.n_val, args.start_init,
		 args.start_comp, args.length_init, args.end, args.k_write, args.primes,args.nmod_primes,args.filename_template_temp,num_thread);
}


void compute_rec_values(fmpz_poly_mat_t Rec,slong r_lcoeff,fmpz_mat_t * Num_inhom,fmpz* Den_inhom,fmpz_mat_t * Num_inits,fmpz* Den_inits,
		bool have_inhom,slong next_prime,slong CHUNK,slong PACKED_CHUNK,slong BITS, slong r_ord,slong rec_degree,slong n_val,slong start_init,
		slong start_comp,slong length_init,slong end,slong k_write,nn_srcptr primes,nmod_t * nmod_primes,const char* filename_template_temp,slong num_thread){
	//nmod_mat_t Values_mod[CHUNK];
	nmod_mat_t * Values_mod=flint_malloc(CHUNK*sizeof(nmod_mat_t));
	nmod_mat_t Values_mod_win,Rec_mod_mat;
	nmod_poly_mat_t Rec_mod;
	nn_ptr k_pow=_nmod_vec_init(rec_degree);
	nn_ptr Rec_mod_k=_nmod_vec_init(r_ord+1);
	nn_ptr Next_val_mod=_nmod_vec_init(n_val);
	nn_ptr num_inhom_k_mod=_nmod_vec_init(n_val);
	nn_ptr val_chunk_k=_nmod_vec_init(CHUNK);
	nn_ptr val_chunk_packed=_nmod_vec_init(PACKED_CHUNK);
	fmpz_t output; fmpz_init(output);
	char filename_temp[300];
	for(slong p=next_prime+num_thread*CHUNK;p<next_prime+(num_thread+1)*CHUNK;p++){//this loop may be parallized 
        //Set modulus to current prime
        nmod_t mod=nmod_primes[p];
        nmod_mat_init(Values_mod[p%CHUNK],n_val,end-start_init+1,primes[p]);
        nmod_mat_init(Rec_mod_mat,r_ord+1,rec_degree,primes[p]);
        //nmod_mat_set_mod(Rec_mod_mat,primes[p]);
        
        nmod_poly_mat_init(Rec_mod, 1, r_ord+1,primes[p]); 
        const dot_params_t rec_eval_dotpar=_nmod_vec_dot_params(rec_degree,mod);
        
        TIMEIT_ONCE_START
        
        
        //reduce Rec modulo mod
        for(slong i=0;i<=r_ord;i++)
            fmpz_poly_get_nmod_poly(nmod_poly_mat_entry(Rec_mod,0,i),fmpz_poly_mat_entry(Rec,0,i));

		for(slong i=0;i<=r_ord;i++)
			for(slong j=0;j<rec_degree;j++)
				nmod_mat_entry(Rec_mod_mat,i,j)=nmod_poly_get_coeff_ui(nmod_poly_mat_entry(Rec_mod,0,i),j);
				
        //Set Values_mod to Values modulo mod
        nmod_mat_window_init(Values_mod_win,Values_mod[p%CHUNK],0,0,n_val,length_init);
        fmpz_mat_get_nmod_mat(Values_mod_win,*Num_inits);
        for(slong k=0;k<length_init;k++){        
			ulong den_inv_mod=nmod_inv(fmpz_get_nmod(Den_inits+k,mod),mod);
			for(slong l=0;l<n_val;l++){
				nmod_mat_entry(Values_mod_win,l,k)=nmod_mul(nmod_mat_entry(Values_mod_win,l,k),den_inv_mod,mod);
			}
        }
        nmod_mat_window_clear(Values_mod_win);
        
        /*//multipoint eval
        nn_ptr * tree=_nmod_poly_tree_alloc(end-start_init+1);
        nn_ptr range=_nmod_vec_init(end-start_init+1);
        nmod_mat_t Evals;
        nmod_mat_init(Evals,r_ord+1,end-start_init+1,mod.n);
        for(slong k=start_init;k<=end;k++)
			range[k-start_init]=k;//reduced!!!!
        _nmod_poly_tree_build(tree,range,end-start_init+1,mod);
        for(slong i=0;i<=r_ord;i++)
			_nmod_poly_evaluate_nmod_vec_fast_precomp(&nmod_mat_entry(Evals,i,0),nmod_poly_mat_entry(Rec_mod,0,i)->coeffs,rec_degree,tree,end-start_init+1,mod);
        _nmod_poly_tree_free(tree,end-start_init+1);
        nmod_mat_clear(Evals);
        _nmod_vec_clear(range);*/
        //------------------------------------
        //             Rec loop
        //------------------------------------
        TIMEIT_ONCE_START
        k_pow[0]=1;
        for(slong k = start_comp; k <= end; k++){
            for(slong i=1;i<rec_degree;i++)
                k_pow[i]=nmod_mul(k_pow[i-1],nmod_sub(k,r_lcoeff,mod),mod);
            slong n_previous=FLINT_MIN(r_ord,k-start_init);
            for(slong i=r_ord-n_previous;i<=r_ord;i++){ //Rec_mod_k=Rec_mod/.n->k
				Rec_mod_k[i]=_nmod_vec_dot(&nmod_mat_entry(Rec_mod_mat,i,0),k_pow,rec_degree,mod,rec_eval_dotpar);
                //Rec_mod_k[i]=poly_evaluate_nmod_powprecomp(nmod_poly_mat_entry(Rec_mod,0,i),k_pow,rec_eval_dotpar);
                //Rec_mod_k[i]=nmod_poly_evaluate_nmod(nmod_poly_mat_entry(Rec_mod,0,i),k); 
            }
            nmod_mat_window_init(Values_mod_win,Values_mod[p%CHUNK],0,k-start_init-n_previous,n_val,k-start_init);
            nmod_mat_mul_nmod_vec(Next_val_mod,Values_mod_win,Rec_mod_k+r_ord-n_previous,n_previous);
            nmod_mat_window_clear(Values_mod_win);
            if(have_inhom){
                for(slong j=0;j<n_val;j++)
				    num_inhom_k_mod[j]=fmpz_get_nmod(fmpz_mat_entry(*Num_inhom,j,k-start_comp),mod); 
				 ulong inverse= nmod_inv(fmpz_get_nmod(Den_inhom+k-start_comp,mod),mod);
				 _nmod_vec_scalar_mul_nmod(num_inhom_k_mod,num_inhom_k_mod,n_val,inverse,mod);
				 _nmod_vec_add(Next_val_mod,Next_val_mod,num_inhom_k_mod,n_val,mod);
			}
                        
            ulong inverse=nmod_inv(Rec_mod_k[r_ord],mod);
            for(slong j=0;j<n_val;j++){
                nmod_mat_entry(Values_mod[p%CHUNK],j,k-start_init)=nmod_mul(Next_val_mod[j],inverse,mod);
            } 
            
        } 
        nmod_mat_clear(Rec_mod_mat);
        nmod_poly_mat_clear(Rec_mod); 
        printf("  loop rec: ");
        TIMEIT_ONCE_STOP;
        //write computed residues to file. Values_mod will be overwritten later to save memory
        if((p+1)%CHUNK==0){
            TIMEIT_ONCE_START               
            FILE * file_tmp= fopen(get_temp_file_name(filename_temp,filename_template_temp,p/CHUNK), "w");
            for(slong k = k_write; k <= end; k++){
                 //flint_printf("write k = %wd\n",k);
                for(slong j=0;j<n_val;j++){ 
                    for(slong l=0;l<CHUNK;l++){
                        val_chunk_k[l]=nmod_mat_get_entry(Values_mod[l],j,k-start_init);
                    }                    
                    _nmod_poly_bit_pack(val_chunk_packed, val_chunk_k, CHUNK, BITS);
                    fmpz_set_ui_array(output,val_chunk_packed,PACKED_CHUNK);
                    //fmpz_set_ui_array(output,val_chunk_k,CHUNK);
                    if(fmpz_out_raw(file_tmp,output)==0){
                        flint_printf("unable to write to file\n");
                    }
                }       
            }
          
            
            fclose(file_tmp);
            printf("  write to file: \n");
            TIMEIT_ONCE_STOP;
        }
        
        printf("computed remainder %lu: ",p);
        TIMEIT_ONCE_STOP;
        //if(end>3260-start_init && n_val>1)
		//	flint_printf("Value: %wu \n",nmod_mat_get_entry(Values_mod[p%CHUNK],1,3260-start_init));
		
    }
    
    for(slong l=0; l<CHUNK;l++)
		nmod_mat_clear(Values_mod[l]);    
    flint_free(Values_mod);
	_nmod_vec_clear(k_pow);
	_nmod_vec_clear(Rec_mod_k);
	_nmod_vec_clear(Next_val_mod);
	_nmod_vec_clear(num_inhom_k_mod);
	_nmod_vec_clear(val_chunk_k);
	_nmod_vec_clear(val_chunk_packed);
	fmpz_clear(output);

}









