#define READ_LEN 250
#define BUFFER_SIZE 1000
#define LINE_LEN 1024
#define GENOME_SIZE 999999
#define FindMax(A,B)    ((A)>(B)?(A):(B))

int GetGenomeLen(char *String){ 
	int x=0,y=0,flag=0, genome_size=0; 
	char temp_loc[20]; 
	for(x=0;x<strlen(String);x++){ 
		if(String[x]=='L' && String[x+1]=='N' && String[x+2]==':'){ 
			x+=3; flag=1; 
		} 
		if(String[x]==' '||String[x]==':'||String[x]=='\n') flag=0; 
		if(flag==1) temp_loc[y++]=String[x]; 
	} 
	temp_loc[y]='\0'; 
	genome_size=atoi(temp_loc); 
	return genome_size;
}

int getLocation(char *String){ 

	int i=0, j=0, tab=0, Location=0;
	char temp_loc[128];

	for(i=0;i<strlen(String);i++){ 
		if(String[i]=='\t') tab++;
		// Take location 
		if(tab==3 && String[i]!='\t') temp_loc[j++]=String[i];
		if(tab==4) {temp_loc[j]='\0'; Location=atoi(temp_loc); break;}
	}
	return Location;
}

int GetInt(char *String, int Column){ 
    int i, j=0, flag=1, IntValue=0; char temp[50]; temp[0]='\0'; 
    for(i=0;;i++){ 
        if(String[i]=='\t' || String[i]=='\n') flag++; 
        if(flag > Column) { temp[j]='\0'; IntValue=atoi(temp);} 
        if(flag==Column && String[i]!='\t' && String[i]!='\n') temp[j++]=String[i]; 
        if(i==strlen(String) || flag > Column) break; 
    } return IntValue;
}

int getStr(char *CharVar, char *String, int Column){ 
	int i, j=0, tab=1; 
	for(i=0;;i++){ 
		if(String[i]=='\t' || String[i]=='\n') tab++; 
		if(tab > Column) { CharVar[j]='\0';} 
		if(tab==Column && String[i]!='\t' && String[i]!='\n') CharVar[j++]=String[i]; 
		if(i==strlen(String) || tab > Column) break; 
	} 
	return strlen(CharVar);
}

void ReWrite(char *CIGAR, char *Seq, char *NewSeq, char *Qual, char *NewQual){ 
	int i=0, j=0, k=0, l=0, m=0, cigar_digit[128]; 
	char cigar_alpha[128], temp_cigar[100];

	// Decode CIGAR 
	for(i=0; i<strlen(CIGAR); i++){

		if(isalpha(CIGAR[i])){ 
			temp_cigar[k]='\0'; 
			cigar_digit[l]=atoi(temp_cigar); 
			cigar_alpha[l]=CIGAR[i]; 
			temp_cigar[0]='\0'; l++; k=0; 
		} 
		if(isdigit(CIGAR[i])){ temp_cigar[k++]=CIGAR[i]; } 
	} 
 	// End of for loop 

	i=j=k=m=0; 

	// rewrite sequence as per cigar 
	for(i=0;i<l;i++){ 
		// soft clip 
		if(cigar_alpha[i]=='S') j+=cigar_digit[i]; 
		// hard clip 
		else if(cigar_alpha[i]=='H') continue; 
		// insertion 
		else if(cigar_alpha[i]=='I') j+=cigar_digit[i]; 
		// deletion 
		else if(cigar_alpha[i]=='D') for(k=j; k<(cigar_digit[i]+j);k++) {NewSeq[m]='-'; NewQual[m]=' '; m++;} 
		// match 
		else if(cigar_alpha[i]=='M'){ 
			for(k=j; k<(cigar_digit[i]+j);k++){
				NewSeq[m]=Seq[k]; 
				NewQual[m]=Qual[k];
				m++; 
			}
		j+=cigar_digit[i]; 
		}
	}
	NewSeq[m]='\0'; NewQual[m]='\0'; // end of rewriting 
}
char GetConsensus(int A, int T, int G, int C, int Aq, int Tq, int Gq, int Cq, int depth, int cutoff){
	char consensus='N';
	int max, Qual=0,total=0;

	total=A+T+G+C;

	if (A+T+G+C >= depth){		// If depth of coverage is equal to or more than asked for
		max=FindMax(FindMax(A,T), FindMax(G,C));
		if(max==A && (A*100)/total >= cutoff) {consensus='A'; Qual=Aq;}
		if(max==T && (T*100)/total >= cutoff) {if(Tq > Qual) {consensus='T'; Qual=Tq;}}
		if(max==G && (G*100)/total >= cutoff) {if(Gq > Qual) {consensus='G'; Qual=Gq;}}
		if(max==C && (C*100)/total >= cutoff) {if(Cq > Qual) {consensus='C'; Qual=Cq;}}
	}
	return consensus;
}
void GetConsensusDetailed(int A,int T,int G,int C,int N,int Aq,int Tq,int Gq,int Cq,int depth,char *consensus,int *coverage,float *af,float *tf,float *gf,float *cf,float *nf,int *aq,int *tq,int *gq,int *cq){
	int max, Qual=0;
	*consensus='N';

	*af=*tf=*gf=*cf=*nf=*aq=*tq=*gq=*cq=0;

	if (A+T+G+C >= depth){		// If depth of coverage is equal to or more than asked for
		max=FindMax(FindMax(A,T), FindMax(G,C));
		if(max==A) {*consensus='A'; Qual=Aq;}
		if(max==T) {if(Tq > Qual) {*consensus='T'; Qual=Tq;}}
		if(max==G) {if(Gq > Qual) {*consensus='G'; Qual=Gq;}}
		if(max==C) {if(Cq > Qual) {*consensus='C'; Qual=Cq;}}
	}

	*coverage=A+T+G+C;

	if(*coverage > 0 ){ 
		if(A>0) {*aq=Aq/A; *aq-=33; *af=A*100/ *coverage;}
		if(T>0) {*tq=Tq/T; *tq-=33; *tf=T*100/ *coverage;} 
		if(G>0) {*gq=Gq/G; *gq-=33; *gf=G*100/ *coverage;} 
		if(C>0) {*cq=Cq/C; *cq-=33; *cf=C*100/ *coverage;}
		if(N>0) {*nf=N*100/ *coverage;}
	} 
}

// This function is used in SAM_SHANNON program
void GetEntropy(int A,int T,int G,int C,int N,int depth,float *entropy){ 
	// Shannon Entropy formula
	
	// H = - SumOf p(x) * log2 p(x);

	int coverage=0; 
	*entropy=0; 
	coverage=A+T+G+C; 
	
	if(coverage > depth ){ 
		if(A>0) *entropy+=(log2((float)A/(float)coverage)*((float)A/(float)coverage)); 
		if(T>0) *entropy+=(log2((float)T/(float)coverage)*((float)T/(float)coverage)); 
		if(G>0) *entropy+=(log2((float)G/(float)coverage)*((float)G/(float)coverage)); 
		if(C>0) *entropy+=(log2((float)C/(float)coverage)*((float)C/(float)coverage)); 
		if(N>0) *entropy+=(log2((float)N/(float)coverage)*((float)N/(float)coverage)); 
		if(*entropy!=0) *entropy=0-*entropy; 
	}
}
