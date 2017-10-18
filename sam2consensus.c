/*
	Author: Sreenu Vattipally
	University of Glasgow, Glasgow

	23/May/2013

	If there are two are more nucleotides with same frequencey, one with the highest qualitied will be picked
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include "sam2consensus.h"

void Usage();

int main(int argc, char **argv){ 

    int i, Location, SeqLen=0, QualLen=0, GenomeLen=0, Flag=0, *A,*T,*G,*C,*N,*Aq,*Tq,*Gq,*Cq,coverage=0, real_genome_size=0;
	char buffer[BUFFER_SIZE], Seq[READ_LEN], NewSeq[READ_LEN], Qual[READ_LEN], NewQual[READ_LEN], CIGAR[READ_LEN], consensus;
	float entropy=0, af=0, tf=0,gf=0,cf=0,nf=0; 
	int aq=0, tq=0,gq=0,cq=0;

	// Command line arguments.
	int q=1, d=1, e=0, out=0,cutoff=0; 
	char InFileName[128], OutFileName[128];

    FILE *input, *output; 
	// Nucleotide Array
	A=(int*) calloc(GENOME_SIZE,sizeof(int));
	T=(int*) calloc(GENOME_SIZE,sizeof(int));
	G=(int*) calloc(GENOME_SIZE,sizeof(int));
	C=(int*) calloc(GENOME_SIZE,sizeof(int));
	N=(int*) calloc(GENOME_SIZE,sizeof(int));
	// Quality Array
	Aq=(int*) calloc(GENOME_SIZE,sizeof(int));
	Tq=(int*) calloc(GENOME_SIZE,sizeof(int));
	Gq=(int*) calloc(GENOME_SIZE,sizeof(int));
	Cq=(int*) calloc(GENOME_SIZE,sizeof(int));

	if(argc < 3){ 
		free(A); free(T); free(G); free(C); free(Aq); free(Tq); free(Gq); free(Cq); free(N);
		Usage();
	}

	// Get command line arguments
	for(i=1; i<argc; i++){
		if(strcmp(argv[i],"-h")==0){ free(A); free(T); free(G); free(C); free(Aq); free(Tq); free(Gq); free(Cq); free(N); Usage();}	// Quality cut-off
		if(strcmp(argv[i],"-q")==0){ q=atoi(argv[i+1]); i++;}	// Quality cut-off
		if(strcmp(argv[i],"-d")==0){ d=atoi(argv[i+1]); i++;}	// Depth cut-off
		if(strcmp(argv[i],"-c")==0){ cutoff=atoi(argv[i+1]); i++;}	// Consensus cut-off
		if(strcmp(argv[i],"-e")==0){ e=1;}					// Extended output
		if(strcmp(argv[i],"-i")==0){ strcpy(InFileName,argv[i+1]); i++;}	// Input file
		if(strcmp(argv[i],"-o")==0){ strcpy(OutFileName,argv[i+1]); i++; out=1;}	// Output file
	}

    if((input = fopen(InFileName, "r"))!=NULL){ 

		while (!feof(input)) { 
			fgets(buffer,LINE_LEN,input); if(feof(input)) break;

				// Get Genome Length
			if(buffer[0]=='@' && strstr(buffer,"LN:")) GenomeLen=GetGenomeLen(buffer);

			if(buffer[0]!='@'){
                Flag=GetInt(buffer, 2); 
                // Take only mapped reads
                if((Flag & 0x0004)!=0x0004){

				Seq[0]='\0'; Qual[0]='\0'; CIGAR[0]='\0'; NewSeq[0]='\0'; NewQual[0]='\0';
				Location=getLocation(buffer);
				getStr(CIGAR, buffer,6);
				SeqLen=getStr(Seq, buffer,10);
				QualLen=getStr(Qual, buffer,11);

				if (Location > GENOME_SIZE-SeqLen){
					printf("Genome size is greater than allocated memory size.\n Please re-compile increasing memory size\n");
					free(A); free(T); free(G); free(C); free(N);
					free(Aq); free(Tq); free(Gq); free(Cq);
					exit(0);
				}

				if(SeqLen==QualLen)
					ReWrite(CIGAR, Seq, NewSeq, Qual, NewQual); 

				//if(Location+strlen(NewSeq) > GenomeLen) GenomeLen=Location+strlen(NewSeq);

				for(i=0; i<strlen(NewSeq); i++){
					if(NewQual[i] >= q+33){					// Consider nucleotides with qualities scores that are equal to or more than asked for
					     if(NewSeq[i]=='A') { A[Location+i-1]++; Aq[Location+i-1]+=NewQual[i]; }
					else if(NewSeq[i]=='T') { T[Location+i-1]++; Tq[Location+i-1]+=NewQual[i]; }
					else if(NewSeq[i]=='G') { G[Location+i-1]++; Gq[Location+i-1]+=NewQual[i]; }
					else if(NewSeq[i]=='C') { C[Location+i-1]++; Cq[Location+i-1]+=NewQual[i]; }
					else N[Location+i-1]++;
					}
				}
			}
        }
		} fclose(input);

		if(out==1) output=fopen(OutFileName,"w"); 

		if(out==1) fprintf(output,">%s\n",InFileName); 
		else printf(">%s\n",InFileName);

	// Printing consensus nucleotides

	for(i=0; i<GenomeLen; i++){
		consensus=GetConsensus(A[i],T[i],G[i],C[i],Aq[i],Tq[i],Gq[i],Cq[i],d,cutoff);
		if(out==1) fprintf(output,"%c",consensus);
		else printf("%c",consensus);
		if(i%70==69){
			if(out==1) fprintf(output,"\n");
			else printf("\n");
		}
	}

	if(out==1) fprintf(output,"\n");
	else printf("\n");

	if(e==1){ 	// Detailed output
		if(out==1) fprintf(output,"Position\tConsensus\tA Freq\tT Freq\tG Freq\tC Freq\tN Freq\tA Qual\tT Qual\tG Qual\tC Qual\tCoveragy\t\tEntropy\n"); 
		else printf("Position\tConsensus\tA Freq\tT Freq\tG Freq\tC Freq\tN Freq\tA Qual\tT Qual\tG Qual\tC Qual\tCoveragy\t\tEntropy\n");

		// Printing positional frequencies of each site
		for(i=0; i<GenomeLen; i++){ 
			GetConsensusDetailed(A[i],T[i],G[i],C[i],N[i],Aq[i],Tq[i],Gq[i],Cq[i], d, &consensus,&coverage,&af,&tf,&gf,&cf,&nf,&aq,&tq,&gq,&cq);
			GetEntropy(A[i],T[i],G[i],C[i],N[i],d,&entropy);
			if (out==1) fprintf (output,"%d\t%c\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t\t%d\t%d\t%d\t%d\t\t%d\t\t%.2f\n", i+1, consensus, af, tf, gf, cf, nf, aq, tq, gq, cq,coverage,entropy); 
			else printf ("%d\t%c\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t\t%d\t%d\t%d\t%d\t\t%d\t\t%.2f\n", i+1, consensus, af, tf, gf, cf, nf, aq, tq, gq, cq,coverage,entropy);
		}
	}
	} free(A); free(T); free(G); free(C); free(Aq); free(Tq); free(Gq); free(Cq); free(N); if(out==1) fclose(output);
}

void Usage(){
	printf("\n\n");
	printf("Program to print consensus nucleotide sequence from a sam file.\n");
	printf("Usage:\n");
	printf("SAM2CONSENSUS -i File.sam\n\n");
	printf("Optional flags\n");
	printf("\t-h This help\n");
	printf("\t-o Output (Output file name. Default: STDOUT)\n");
	printf("\t-d integer (Minimum depth of coverage. default: 1)\n");
	printf("\t-q integer (Minimum PHRED quality of nucleotide.  default: 1)\n");
	printf("\t-e (Extended output  default: NO)\n");
	printf("\n\n");
	exit (0);
}
