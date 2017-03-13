#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define MAXNUM 200
#define MAXLINE 600000

struct gene{
  char filn;
  float chr;
  int sta;
  int ter;
  int lin;
};


void freegenes (struct gene * l[], int i){
  int j;
  for(j=0; j<i; j++)
    free(l[j]);
}

void strcopy(char s[],char t[]);

/* readfile: input DNA fragment position data */
int readfile(FILE * fp, int startline, struct gene * l[], int s, int e){
  int i=0, j=0, r, k, d=0, n=0;
  char t[50];
  char w,v; 
  w=fgetc(fp);
  while ( startline > 1 ){
    if((v=fgetc(fp))=='\n' && w!='\n')
      startline--;
    w=v;
    }

  l[j]=(struct gene * )malloc(sizeof(struct gene));
  l[j]->chr=0;
  while ( (k=fgetc(fp))!=EOF ){
    if ( k=='\n' && i!=0 ){
      l[j]->lin=j;
      j++;
      if (j==MAXLINE)
        break;
      l[j]=(struct gene * )malloc(sizeof(struct gene));
      l[j]->chr=0;
      i=0;
      n=0;
      t[n]='\0';
    }
    else if ( k=='+' ){
      l[j]->chr += 0;
      t[n++]=k;}
    else if( k==' ' || k=='\t' || k==','){
      if(n==0)
        continue;
      i++;
      t[n]='\0';
      n=0;
      if ( t[0]=='c'&&t[1]=='h'&&t[2]=='r' ){
        for(r=3;r<strlen(t);r++){
          d=d*10+(t[r]-'0');
        }
        l[j]->chr += d;
        d = 0;
      }
      else if (i==s)
        l[j]->sta=atoi(t);
      else if (i==e)
        l[j]->ter=atoi(t);
    }
    else{
      t[n++]=k;
      
    }
  }
  return j;
}

/* strcopy: copy string s to string t */
void strcopy(char s[],char t[]){
  int i,j;
  for (i=0,j=0;j<strlen(t);j++){
    if(t[j]=='"' || t[j]==';')
      ;
    else
      s[i++]=t[j];
  }
  s[i]='\0';

}


/* swap: interchange v[i] and v[j] */ 
void swap(struct gene *v[], int i, int j) 
{
  struct gene *temp; 

  temp = v[i]; 
  v[i] = v[j]; 
  v[j] = temp; 
}




void qsortposi (struct gene  *v[], int left, int right) {
  int i, last; 
  void swap (struct gene *v[], int i, int j); 
  if (left >= right) /* do nothing if array contains */
    return;          /* fewer than two elements */ 
  swap(v, left, (left + right)/2);  /* move partition element to v[0] */ 
  last = left; 
  for ( i = left + 1; i <= right; i++) /* partition */ 
    if ((v[i]->sta < v[left]->sta && v[i]->chr == v[left]->chr) || v[i]->chr < v[left]->chr     ) 
      swap(v, ++last, i); 
  swap(v, left, last); /* restore partition elem */
  qsortposi ( v, left, last - 1); 
  qsortposi (v, last +1, right); 
}


void qsortposie (struct gene  *v[], int left, int right) {
  int i, last; 
  void swap (struct gene *v[], int i, int j); 
  if (left >= right) /* do nothing if array contains */
    return;          /* fewer than two elements */ 
  swap(v, left, (left + right)/2);  /* move partition element to v[0] */ 
  last = left; 
  for ( i = left + 1; i <= right; i++) /* partition */ 
    if ( (v[i]->ter < v[left]->ter &&  v[i]->chr == v[left]->chr) || v[i]->chr < v[left]->chr) 
      swap(v, ++last, i); 
  swap(v, left, last); /* restore partition elem */
  qsortposie ( v, left, last - 1); 
  qsortposie (v, last +1, right); 
}
 

void qsortlin (struct gene  *v[], int left, int right) {
  int i, last;
  void swap (struct gene *v[], int i, int j);
  if (left >= right) /* do nothing if array contains */
    return;          /* fewer than two elements */
  swap(v, left, (left + right)/2);  /* move partition element to v[0] */
  last = left;
  for ( i = left + 1; i <= right; i++) /* partition */
    if (v[i]->lin < v[left]->lin )
      swap(v, ++last, i);
  swap(v, left, last); /* restore partition elem */
  qsortlin ( v, left, last - 1);
  qsortlin (v, last +1, right);
}



/* compare two DNA fragments position and decide overlapping DNA fragments */
void comgeneposi (struct gene *la[], int numa, struct gene * lb[], int numb){
  int i=0,j=0,star=0;
  float k=0.0;
  for (i=0; i<numa ; i++){
    j = star;
    if(la[i]->chr < lb[j]->chr || (la[i]->chr == lb[j]->chr && la[i]->ter < lb[j]->sta)){
      if(la[i]->filn != 'c'){
	la[i]->filn = 'a';
      }
      continue;
    }
    while ((la[i]->sta > lb[j]->ter && la[i]->chr == lb[j]->chr) || la[i]->chr > lb[j]->chr){
      if(lb[j]->filn!='d')
        lb[j]->filn='b';
      j++;
      if (j==numb)
        break;
    }
    if (j==numb)
      break;
    if (j==0){
      star = j;
    }else{
      star = j - 1;
    }
    
    while (la[i]->ter >= lb[j]->sta && la[i]->chr==lb[j]->chr){
		if(lb[j]->sta <= la[i]->sta && la[i]->ter <= lb[j]->ter ){
			lb[j]->filn='d';
			k=1.0;
		}else if(lb[j]->sta <= la[i]->sta && la[i]->sta <= lb[j]->ter){
			lb[j]->filn='d';
			k = k + 1.0*(lb[j]->ter - la[i]->sta)/(la[i]->ter - la[i]->sta);
		}else if(la[i]->sta < lb[j]->sta && la[i]->ter < lb[j]->ter){
			lb[j]->filn='d';
			k = k + 1.0*(la[i]->ter - lb[j]->sta)/(la[i]->ter - la[i]->sta);
		}else if(la[i]->sta < lb[j]->sta && la[i]->ter >= lb[j]->ter){
			lb[j]->filn='d';
			k = k + 1.0*(lb[j]->ter - lb[j]->sta)/(la[i]->ter - la[i]->sta);
		}else{
			lb[j]->filn='b';
		}
        j++;
        if(j==numb)
          break;
    }
    if ( k > 0.1 )
      la[i]->filn='c';
    else 
      la[i]->filn='a';
    k=0.0;
  }
  for( ; i<numa; i++)
    la[i]->filn='a';
  for( ; j<numb; j++)
    if(lb[j]->filn!='d')
      lb[j]->filn='b';
}

/* output overlapping DNA fragments */
int pfile (char filea[], int j, struct gene * l[], int i){
  FILE * fa, * fb, * fab_a, * fab_b, * fp1;
  char t[500];
  fa=fopen("A-B","a");
  fb=fopen("B-A","a");
  fab_a=fopen("A&B_A","a");
  fab_b=fopen("A&B_B","a");
  fp1=fopen(filea,"r");
  int r, k=0;
  while ( ++k < j ){
    fgets(t,500,fp1);
    while(t[0]=='\n')
      fgets(t,500,fp1);
    }
  for (r=0; r<i; r++){
    fgets(t,500,fp1);
    while ( t[0]=='\n' )
      fgets(t,500,fp1);
    if(l[r]->filn=='a')
      fputs(t,fa);
    else if(l[r]->filn=='b')
      fputs(t,fb);
    else if(l[r]->filn=='c')
      fputs(t,fab_a);
    else
      fputs(t,fab_b);
  }
  fclose(fa);
  fclose(fb);
  fclose(fab_a);
  fclose(fab_b);
  fclose(fp1);
  return i+j;
}


int recofile (int j, struct gene * l[], int i){
  FILE * fre, * tem;
  fre = fopen("record_file.txt","r");
  tem = fopen("temple.txt","w");
  char t[500];
  int k = 1, len;
  char ch;
  while (k < j){
    fgets(t, 500, fre);
    fputs(t, tem);
        k++;
  }
  k = 0;
  while (k < i){
    fgets(t, 500, fre);
    len = strlen(t);
    if(l[k]->filn=='a'){
      t[len-1]='a';
      t[len]='\n';
      t[len+1]='\0';
    }
    else{
      t[len-1]='n';
      t[len]='\n';
      t[len+1]='\0';
    }
    fputs(t, tem);
        k++;
  }
  k = i+j;
  while (k <= 200000){
    fgets(t, 500, fre);
    fputs(t, tem);
    k++;
  }
  fclose(tem);
  fclose(fre);
  fre = fopen("record_file.txt","w");
  tem = fopen("temple.txt","r");
  ch = fgetc(tem);
  while(ch != EOF){
    fputc(ch, fre);
    ch = fgetc(tem);
  }
  fclose(tem);
  fclose(fre);
}


int main(int argc, char *argv[]){
  clock_t startt, endt;
  float tim;
  startt = clock(); 
  FILE *fp1, *fp2;
  struct gene * la[MAXLINE], * lb[MAXLINE];
  int i, j, t=0, k=1, s=1, tem;
  int starta=0, startb=0, enda=0, endb=0;
  if (argc !=8 ){
    printf("usage:Biodiff -c -a 3,4 -b 3,4 A_ucsc_genes.txt B_ucsc_gene.gtf\n");
    return 0;
  }
  if (strcmp(argv[1],"-c")==0){
    while(argv[3][t]!=',')                  
      starta=starta*10 + argv[3][t++]-'0'; 
    while(argv[3][++t]!='\0')
      enda=enda*10 + argv[3][t]-'0';
    t=0;
    while(argv[5][t]!=',')
      startb=startb*10 + argv[5][t++]-'0';
    while(argv[5][++t]!='\0')
      endb=endb*10 + argv[5][t]-'0';
    if((fp1=fopen(argv[6],"r"))==NULL){
      printf("Can not find the %s\n",argv[6]);
      return 0;}
    if((fp2=fopen(argv[7],"r"))==NULL){
      printf("Can not find the %s\n",argv[7]);
      return 0;}
    while( (i = readfile(fp1,s,la,starta,enda))==MAXLINE ){
                fclose(fp1);
                qsortposi(la,0,i-1);
                k=1;
                while( (j=readfile(fp2,k,lb,startb,endb))==MAXLINE ){
                fclose(fp2);
                qsortposie(lb,0,j-1);
                comgeneposi(la,i,lb,j);
                qsortlin(lb,0,j-1);
                k=pfile(argv[7],k,lb,j);
                fp2=fopen(argv[7],"r");
                freegenes(lb,j);
                printf("%d\t%d\n",k,j);
                }
                fclose(fp2);
   
                qsortposie(lb,0,j-1);
                comgeneposi(la,i,lb,j);
                qsortlin(la,0,i-1);
                qsortlin(lb,0,j-1);
                recofile(s,la,i);
                printf("record ok!\n");
                s = pfile(argv[6],s,la,i);
                k = pfile(argv[7],k,lb,j);
                fp1=fopen(argv[6],"r");
                fp2=fopen(argv[7],"r");
                freegenes(la,i);
                freegenes(lb,j);
                printf("%d\n>%d\n",k,s);
        }
        fclose(fp1);
        qsortposi(la,0,i-1);
        k=1;
        while( (j=readfile(fp2,k,lb,startb,endb))==MAXLINE ){
                fclose(fp2);
                qsortposie(lb,0,j-1);
                comgeneposi(la,i,lb,j);
                qsortlin(lb,0,j-1);
                k=pfile(argv[7],k,lb,j);
                fp2=fopen(argv[7],"r");
                freegenes(lb,j);
                printf("%d\n",k);
        }
        fclose(fp2);
   
        qsortposie(lb,0,j-1);
        comgeneposi(la,i,lb,j);
        qsortlin(la,0,i-1);
        qsortlin(lb,0,j-1);
        recofile(s,la,i);
        s = pfile(argv[6],s,la,i);
        k = pfile(argv[7],k,lb,j);
        printf("%d\n>%d\n",k,s);
        freegenes(la,i);
        freegenes(lb,j);
  }

  endt = clock();
  tim = (endt-startt)/CLOCKS_PER_SEC;
  printf("compare gene between %s and %s cost %f second.\n", argv[6],argv[7],tim);
  return 1;
}
