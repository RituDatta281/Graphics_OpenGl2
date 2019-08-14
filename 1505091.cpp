#include<iostream>
#include<cstring>
#include<cstdlib>
#include<fstream>
#include<stack>
#include<cmath>
#include <iomanip> 
#define pi 3.14159265 
using namespace std;
ofstream output1("stage1.txt");
ofstream output2("stage2.txt");
ofstream output3("stage3.txt");
struct point
{
	double x,y,z;
}scale,translate, triangle[3];
struct rotation
{
	double angle,x,y,z;
}rotate;

struct Matrix
{
	double matrix[4][4];
	Matrix()
	{
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				matrix[i][j]=0;
		matrix[3][3]=1;
	}

} identity, Current,V_view,projectionMat;

double fovY,aspectRatio, near,far;
string function;
struct vector_
{
	double I,J,K;
	vector_(double x,double y, double z)
	{
		I=x;
		J=y;
		K=z;
	}
	vector_()
	{

	}

}eye,look,up;

stack<Matrix> mystack; 
stack<int> index_;

Matrix normalise(Matrix a)
{
	if(a.matrix[3][3]!=1)
	{
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				a.matrix[i][j]= a.matrix[i][j]/a.matrix[3][3];
	}
	
	return a;
}
void printMatrix( double matrix[4][4])
{
	// matrix= new double*[row];
	// for(int i=0;i<row;i++)
	// 	matrix[i]= new double[column];

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
			{
				cout<<setprecision(7)<<matrix[i][j]<<" ";
			}
			cout<<endl;
	}
}

Matrix Matrix_multiplication(Matrix a, Matrix b)
{
	Matrix result;
	 for (int i = 0; i < 4; i++) 
    { 
        for (int j = 0; j < 4; j++) 
        { 
            result.matrix[i][j] = 0; 
            for (int k = 0; k < 4; k++) 
                result.matrix[i][j] += a.matrix[i][k] *  
                             b.matrix[k][j]; 
        } 
    } 
    return result;
}

double dotProduct(vector_ a, vector_ b)
{
	return a.I*b.I+ a.J*b.J+ a.K*b.K;
}

vector_ crossProduct(vector_ u, vector_ v)
{
	vector_ result;
	result.I=u.J* v.K - v.J * u.K;
	result.J=v.I * u.K - u.I * v.K;
	result.K=u.I * v.J - v.I * u.J;
	return result;
}



void multiply(double a[][4],double b[][1], ofstream &f)
{
	 double res[4][1]; 

	 for (int i = 0; i < 4; i++) 
    { 
        for (int j = 0; j < 1; j++) 
        { 
            res[i][j] = 0; 
            for (int k = 0; k < 4; k++) 
                res[i][j] += a[i][k] * b[k][j]; 
        } 
    } 

    for (int i = 0; i < 1; i++)  
    { 
        for (int j = 0; j < 4-1; j++)  
        { 	

           f<< std::fixed << std::setprecision(7)<< res[j][i]/ res[3][i] << " "; 
        } 
       // f << "\n"; 
    } 
    f << "\n"; 
    
}




void vectorPrint(vector_ v)
{
	cout<<v.I<<"-"<<v.J<<"-"<<v.K<<endl;
}

vector_ normalise_vector(vector_ a)
{
	double r = sqrt(a.I*a.I + a.J*a.J + a.K*a.K);
	a.I = a.I/r;
	a.J = a.J/r;
	a.K = a.K/r;
	//vectorPrint(a);
	return a;

}

void init()
{

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
			{
				if(i==j) identity.matrix[i][j]=1;
				else identity.matrix[i][j]=0;
			}
	}


//	printMatrix(identity.matrix);

	mystack.push(identity);

	//printMatrix(mystack.top().matrix);
	



}
vector_ R_(vector_ x,vector_ a,double angle)
{
	vector_ result;
	if(angle==90.0) angle=pi/2;
	else
	angle= (angle*pi)/180.0;
	float val=cos(angle);
	result.I= cos(angle)*x.I+ (1-cos(angle))*dotProduct(a,x)*a.I+ (sin(angle)*crossProduct(a,x).I);
	result.J= cos(angle)*x.J+ (1-cos(angle))*dotProduct(a,x)*a.J+ (sin(angle)*crossProduct(a,x).J);
	result.K= cos(angle)*x.K+ (1-cos(angle))*dotProduct(a,x)*a.K+ (sin(angle)*crossProduct(a,x).K);

	//cout<<setprecision(8)<<val<<" "<<sin(angle)<<endl;

	//printf("%0.2f\n",cos(angle) );

	//cout<<(sin(angle)*crossProduct(a,x).I)<<endl;

	return result;


}

void projectionTransformation()
{
	double fovX = fovY * aspectRatio;
	double t = near * tan((fovY/2)*pi/180);
	double r = near * tan((fovX/2)*pi/180);

	projectionMat.matrix[0][0]=near/r;
	projectionMat.matrix[1][1]=near/t;
	projectionMat.matrix[2][2]=-(far+near)/(far-near);
	projectionMat.matrix[2][3]=-(2*far*near)/(far-near);
	projectionMat.matrix[3][3]=0;
	projectionMat.matrix[3][2]=-1;


}





void viewTransformation()
{
	vector_ l(look.I-eye.I,look.J-eye.J,look.K-eye.K);
	//vectorPrint(l);
	l=normalise_vector(l);
	//vectorPrint(l);
	vector_ r= crossProduct(l,up);
	//vectorPrint(r);
	normalise_vector(r);
	vector_ u= crossProduct(r,l);

	Matrix T_view,R_view;
	for(int i=0;i<3;i++)
	{

		for(int j=0;j<3;j++)
		{
			if(i==j)
				T_view.matrix[i][j]=1;
		}
	}

	T_view.matrix[0][3]=-eye.I;
	T_view.matrix[1][3]=-eye.J;
	T_view.matrix[2][3]=-eye.K;

	R_view.matrix[0][0]=r.I;
	R_view.matrix[0][1]=r.J;
	R_view.matrix[0][2]=r.K;
	R_view.matrix[1][0]=u.I;
	R_view.matrix[1][1]=u.J;
	R_view.matrix[1][2]=u.K;
	R_view.matrix[2][0]=-l.I;
	R_view.matrix[2][1]=-l.J;
	R_view.matrix[2][2]=-l.K;


	V_view=Matrix_multiplication(R_view,T_view);
	// cout<<"VIEW"<<endl;
	printMatrix(V_view.matrix);

 }


void stage2()
{
	ifstream fin("stage1.txt");
	double f;
	int count=0;
	while(fin>>f)
	{
	
		double triangleMat[4][1];
		triangleMat[0][0]=f;
		fin>>triangleMat[1][0]>>triangleMat[2][0];
		triangleMat[3][0]=1;

	  	multiply(V_view.matrix,triangleMat,output2);

	  	count++;
	   	if(!(count%3))
	   		output2<<endl;
	}
	 output2.close();
	 fin.close();	
}

void stage3()
{
	ifstream fin("stage2.txt");
	double f;
	int count=0;
	while(fin>>f)
	{
	
		double triangleMat[4][1];
		triangleMat[0][0]=f;
		fin>>triangleMat[1][0]>>triangleMat[2][0];
		triangleMat[3][0]=1;

	  	multiply(projectionMat.matrix,triangleMat,output3);

	  	count++;
	   	if(!(count%3))
	   		output3<<endl;
	}
	 
}



int main()
{
	ifstream input;
	input.open("scene.txt");
	
	string line;

	input>>eye.I>>eye.J>>eye.K;
	input>>look.I>>look.J>>look.K;
	input>>up.I>>up.J>>up.K;
	input>>fovY>>aspectRatio>>near>>far;
	//cout<<fovY<<" "<<aspectRatio<<" "<<near<<" "<<far<<endl;


	init();

	while(input>>function)
	{
	
	if(!function.compare("triangle"))
	{
		//point triangle;
		for(int i=0;i<3;i++)
		{
		double triangleMat[4][1];

		input>>triangleMat[0][0]>>triangleMat[1][0]>>triangleMat[2][0];
		triangleMat[3][0]=1;

	  	multiply(mystack.top().matrix,triangleMat, output1);
	  	}
	  	output1<<endl;
	}

	if(!function.compare("push"))
	{
		index_.push(mystack.size());
	}

	if(!function.compare("pop"))
	{
		int index=index_.top();
		while(mystack.size()!=index)
		{
			mystack.pop();
		}

		index_.pop();
	}

	 if(!function.compare("scale"))
	 {
	 	Matrix scalingMat;

	 	for(int i=0;i<3;i++)
	 		for(int j=0;j<3;j++)
	 		{
	 			if(i==j)
	 				input>>scalingMat.matrix[i][j];
	 		}

	 		Matrix previous=mystack.top();
	 		Matrix newMat= Matrix_multiplication(previous,scalingMat);
	 		Matrix new_Mat= normalise(newMat);
	 	//	printMatrix(new_Mat.matrix);
	 		mystack.push(new_Mat);


	 }

	if(!function.compare("translate"))
	{
		Matrix translationMat;

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				if(i==j)
					translationMat.matrix[i][j]=1;

		input>>translationMat.matrix[0][3]>>translationMat.matrix[1][3]>>translationMat.matrix[2][3];
			Matrix previous=mystack.top();
	 		Matrix newMat= Matrix_multiplication(previous,translationMat);
	 		Matrix new_Mat= normalise(newMat);
	 	//	printMatrix(new_Mat.matrix);
	 		mystack.push(new_Mat);
	}

	if(!function.compare("rotate"))
	{
		input>>rotate.angle>>rotate.x>>rotate.y>>rotate.z;
		vector_ i_(1.0,0.0,0.0);
		vector_ j_(0.0,1.0,0.0);
		vector_ k_(0.0,0.0,1.0);

		vector_ a(rotate.x,rotate.y,rotate.z);

		a=normalise_vector(a);

		vector_ C1,C2,C3;

		C1=R_(i_,a,rotate.angle);
		C2=R_(j_,a,rotate.angle);
		C3=R_(k_,a,rotate.angle);

		Matrix RotaionMatrix;

		RotaionMatrix.matrix[0][0]=C1.I;
		RotaionMatrix.matrix[1][0]=C1.J;
		RotaionMatrix.matrix[2][0]=C1.K;
		RotaionMatrix.matrix[0][1]=C2.I;
		RotaionMatrix.matrix[1][1]=C2.J;
		RotaionMatrix.matrix[2][1]=C2.K;
		RotaionMatrix.matrix[0][2]=C3.I;
		RotaionMatrix.matrix[1][2]=C3.J;
		RotaionMatrix.matrix[2][2]=C3.K;

		//printMatrix(RotaionMatrix.matrix);

			Matrix previous=mystack.top();
	 		Matrix newMat= Matrix_multiplication(previous,RotaionMatrix);
	 		Matrix new_Mat= normalise(newMat);
	 		printMatrix(new_Mat.matrix);
	 		mystack.push(new_Mat);

	}

	}


	output1.close();
	input.close();  


	viewTransformation();
	stage2();


	projectionTransformation();
	stage3();



	return 0;
}