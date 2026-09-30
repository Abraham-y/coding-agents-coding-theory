/* complement_edges.c -- for a nearly complete mis_ils graph, emit the NON-edges
 * (the compatibility graph); an independent set of the input is a clique here.
 * usage: complement_edges in.txt out.txt   (out: k m' / weights / non-edges) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc,char**argv){ if(argc<3)return 2; FILE*f=fopen(argv[1],"r"); long k,m; if(fscanf(f,"%ld %ld",&k,&m)!=2)return 2;
 long*w=malloc(sizeof(long)*k); for(long i=0;i<k;i++) if(fscanf(f,"%ld",&w[i])!=1)return 2;
 size_t W=(size_t)k*(size_t)k; uint64_t*b=calloc((W+63)/64,8); long a,c;
 while(fscanf(f,"%ld %ld",&a,&c)==2){ size_t i=(size_t)a*k+c; b[i>>6]|=1ULL<<(i&63); i=(size_t)c*k+a; b[i>>6]|=1ULL<<(i&63);} fclose(f);
 long cnt=0; for(long i=0;i<k;i++)for(long j=i+1;j<k;j++){size_t x=(size_t)i*k+j; if(!((b[x>>6]>>(x&63))&1))cnt++;}
 FILE*o=fopen(argv[2],"w"); fprintf(o,"%ld %ld\n",k,cnt); for(long i=0;i<k;i++)fprintf(o,"%ld%s",w[i],i+1<k?" ":"\n");
 for(long i=0;i<k;i++)for(long j=i+1;j<k;j++){size_t x=(size_t)i*k+j; if(!((b[x>>6]>>(x&63))&1))fprintf(o,"%ld %ld\n",i,j);} fclose(o);
 fprintf(stderr,"%ld vertices, %ld non-edges\n",k,cnt); return 0;}
