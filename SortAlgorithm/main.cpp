#include <iostream>
#include "stdlib.h"
using namespace std;
class Sort_Algo{
public:
    void insertsort(int a[],int n)
    {
        int j;
        for(int i=1;i<n;i++)
        {
            if(a[i]<a[i-1])
            {
                int tmp=a[i];
                for(j=i-1;j>=0&&a[j]>tmp;j--)
                {
                    a[j+1]=a[j];
                }
                a[j+1]=tmp;
            }
        }
    }
     void shellsort(int a[],int n)
    {
        int i,j,tmp;
        for(int d=n/2;d>=1;d=d/2){
            for(i=d;i<n;i++){
                tmp=a[i];
                for(j=i-d;j>=0&&tmp<a[j];j=j-d){
                    a[j+d]=a[j];
                }
                a[j+d]=tmp;
            }
        }
    }
    void selectsort(int a[],int n){
        int tmp=0;
        int k;
        for(int i=0;i<n;i++){
            k=i;
            for(int j=i+1;j<n;j++){
                if(a[k]>a[j]){
                    k=j;
                }
            }
            if(k!=i){
                tmp=a[i];
                a[i]=a[k];
                a[k]=tmp;
            }
        }
    }

    void heapadjust(int a[],int i,int n){
        int nchild;
        int ntemp;
        for(;2*i+1<n;i=nchild){
            nchild=2*i+1;
            if(nchild<n-1&&a[nchild+1]>a[nchild])
                ++nchild;
            if(a[i]<a[nchild]){
                ntemp=a[i];
                a[i]=a[nchild];
                a[nchild]=ntemp;
            }else{
                break;
            }
        }
    }
    void heapsort(int a[],int n){//¶ÑÅÅÐò
        int tmp;
        for(int i=n/2-1;i>=0;--i){
            heapadjust(a,i,n);
            for(int j=n-1;i>0;--j){
                tmp=a[i];
                a[i]=a[0];
                a[0]=tmp;
                heapadjust(a,0,i);
            }
        }
    }




    void bubblesort(int a[],int n){
        int i,j,t;
        for(i=0;i<n-1;i++){
           for(j=0;j<n-1-i;j++){
                if(a[j]>a[j+1]){
                    t=a[j];
                    a[j]=a[j+1];
                    a[j+1]=t;
                }

           }
        }

    }

    void bubblesort1(int a[],int n){
        int i ,j ,t;
        for(i=0;i<n;i++)
            for(j=i+1;j<n;j++)
                if(a[i]>a[j])
                {
                    t=a[i];
                    a[i]=a[j];
                    a[j]=t;
                }
        }


    void quicksort(int a[],int low,int high){
        if(low<high){
            int first=low;
            int last=high;
            int key=a[first];
            while(first<last){
                while(first<last&&a[last]>=key){
                    last--;
                }
                a[first]=a[last];
                while(first<last&&a[first]<=key){
                    first++;
                }
                a[last]=a[first];
            }
            cout<<last<<first<<endl;
            a[first]=key;
            quicksort(a,low,first-1);
            quicksort(a,first+1,high);
        }


    }
};
void test(int low)
{
    low++;

    cout<<"low:"<<low<<endl;
    if(low<10){
        test(low);
    }
    //cout<<"high:"<<high<<endl;

}
int main()
{
    //test(1);
    int a[]={4,9,5,12,12,99,4,78,45,22,1,8,9,66};
    int n=sizeof(a)/sizeof(a[0]);
   // int b[14];
    Sort_Algo sort_method;
    sort_method.bubblesort1(a,n);
    //sort_method.MergeSort(a,b,0,13);
   // is.insertsort(a,n);
    //ssort ss;
   // ss.shellsort(a,n);
    for(int index=0;index<n;index++)

       cout<<a[index]<<" ";

    system("pause");
    return 0;
}
