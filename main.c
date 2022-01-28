#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define NB_OBJ 7
typedef struct t_chromosome {
    unsigned char obj;
    int F,L;
}chromosome;
chromosome* generer_chromosome();
void disp_chromosome(chromosome* C);
int getb(unsigned char octet,int n);
void setb(unsigned char* octet,int n);
void shiftb(unsigned char* octet,int n);
void muter(chromosome *C);
void trier_tab_chromosomes(chromosome** tab,int sz);
void tirer_ind_chr(chromosome** tab,int sz, int* p1, int* p2);
void newGen(chromosome** tab,int p1, int p2, int sz);
void penaliser(chromosome* C);
void suppr_chr(chromosome* C);
void blank();
void tprintf(const char* msg,int insert);
void tdisp_chromosome(chromosome* C);
double U01();
int randint(int a,int b);
int m = 7;
int penal = 1000;
int NB_GEN = 200;
int V[NB_OBJ+1] = {0,70,20,39,37,7,5,10};
int W[NB_OBJ+1] = {0,31,10,20,19,4,3,6};
int b  = 50;
trace = 1;
int main(int argc, char* argv[])
{
    srand(time(NULL));
    chromosome** tab_chr=calloc(m+1,sizeof(chromosome*));
    for(int i=1; i<=m; tab_chr[i]=generer_chromosome(),i++);
    trier_tab_chromosomes(tab_chr,m);
    printf("GENERATION 1 :\n");
    puts("__________________________");
    for(int i=1; i<=m; disp_chromosome(tab_chr[i]),i++);
    for(int i=1; i<=NB_GEN; i++){
        int p1,p2;
        tirer_ind_chr(tab_chr,m,&p1,&p2);
        newGen(tab_chr,p1,p2,m);
        if(trace || i==NB_GEN){
            printf("\nGENERATION %d :\n",i+1);
            puts("__________________________");
            for(int i=1; i<=m; disp_chromosome(tab_chr[i]),i++);
        }
    }
    printf("\nOPTIMUM :\n");
    puts("__________________________");
    disp_chromosome(tab_chr[1]);

	for(int i=1;i<=m; suppr_chr(tab_chr[i]),i++);
	free(tab_chr);
    return 0;
}

double U01(){
    return rand()/(double)RAND_MAX;
}
int randint(int a,int b){
    return (a+(int)floor(U01()*(b-a)));
}

chromosome* generer_chromosome() {
    chromosome* C=calloc(1,sizeof(chromosome));
    int nc = 0;
    int randi=0;
    int tab_dispo[NB_OBJ+1];
    for(int i=0; i<NB_OBJ+1;tab_dispo[i]=i,i++);
    int n;
    do{
        nc++;
        randi = randint(1,NB_OBJ);
        n = 0;
        while((tab_dispo[randi]==0 || (tab_dispo[randi]!=0)&&C->L + W[randi] > b )&& n<NB_OBJ) {
                randi=(randi+1)%NB_OBJ + 1;
                n++;
        }
        if(!(n>=NB_OBJ)){
            setb(&(C->obj),randi);
            C->F+= V[randi];
            C->L+=W[randi];
        }
    }while(n<NB_OBJ);
    return C;
}

void disp_chromosome(chromosome* C){
    for(int i=1; i<=NB_OBJ; printf("%d ",getb(C->obj,i)),i++);
    printf(", F = %d , L = %d\n",C->F,C->L);
}
void tdisp_chromosome(chromosome* C){
    (trace)?disp_chromosome(C):blank();
}
void trier_tab_chromosomes(chromosome** tab,int sz){
    int bestchr=1;
    for(int i =2; i<=sz; i++){
        if(tab[bestchr]->F<tab[i]->F) bestchr=i;
    }
    chromosome* tmp = tab[1];
    if(bestchr!=1){
        tab[1]=tab[bestchr];
        for(int i=bestchr;i>2;i--){
            tab[i]=tab[i-1];
        }
        tab[2]=tmp;
    }
    return;
}

void tirer_ind_chr(chromosome** tab,int sz, int* p1, int* p2){
    int p=randint(1,sz);
    int q;
    do { q = randint(1,sz); } while (q==p);
    if (tab[p]->F >= tab[q]->F) *p1 = p;
    else *p1 = q;
    do {*p2 = randint(1,sz);} while(*p2==*p1);

    return;
}
void suppr_chr(chromosome* C){
    free(C);
	return;
}

void newGen(chromosome** tab,int p1, int p2,int sz){
	chromosome P1 = *tab[p1],P2 = *tab[p2];
	tprintf("->Parent P1 :\n\t",0); tdisp_chromosome(&P1);
	tprintf("->Parent P2 :\n\t",0); tdisp_chromosome(&P2);
	tputchar('\n');
	chromosome* newChr;
	for(int k=2;k <=sz;k++){
		newChr=calloc(1,sizeof(chromosome));
		int cut = 2+(int)floor(U01()*6);
        tprintf("\t->Cut : %d\n",cut);
		for(int i =1 ; i<cut ; i++) {
		    if(getb(P1.obj,i)){
		        setb(&(newChr->obj),i);
		        newChr->F+=V[i];
		        newChr->L+=W[i];
		    }
		}
		for(int i =cut ; i<=NB_OBJ ; i++) {
		    if(getb(P2.obj,i)){
		        setb(&(newChr->obj),i);
		        newChr->F+=V[i];
		        newChr->L+=W[i];
		    }
		}
		muter(newChr);
		penaliser(newChr);
        tprintf("\t->Enfant E%d : ",k-1); tdisp_chromosome(newChr);
        tputchar('\n');
		suppr_chr(tab[k]);
		tab[k]=newChr;
	}
	trier_tab_chromosomes(tab,sz);
	tputchar('\n');
    return;
}

void penaliser(chromosome* C){
    if(C->L>b) {
        C->F-=(C->L-b)*penal;
        tprintf("\t-> penalite : on inflige -%d\n",(C->L-b)*penal);
    }
    return;
}
void muter(chromosome* C){
    int n = randint(1,NB_OBJ);
	shiftb(&(C->obj),n);
	tprintf("\t-> Mutation du bit %d\n",n);
	return;
}

int getb(unsigned char octet,int n){
    return (octet & 1<<n)>>n;
}
void setb(unsigned char* octet,int n){
    *octet=(*octet)|(1<<n);
}
void shiftb(unsigned char* octet,int n){
    *octet+=(getb(*octet,n))?-(1<<n):(1<<n);
	return;
}
void blank(){
    printf("");
    return;
}
void tprintf(const char* msg,int insert){
    char msg2[90];
    sprintf(msg2,msg,insert);
    (trace)?printf("%s",msg2):blank();
    return;
}
void tputchar(char c){
    (trace)?putchar(c):blank();
}
