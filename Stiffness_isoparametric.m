clear all
close all
clc
%%
%% element stiffness matrix for quad element using 4 point integration
%% material data
E=128.5*10^9;%Gpa to pa
v=0.119%0.25;
t=0.210;%m
%% Gauss weights
W1=1
W2=1
int_points=4;

nodes=4
%% gauss points

zi =[-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
eta=[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];

%% Global coordinates       
 x1= -1.9000000;x2=5.7000000;x3=7.6000000;x4=-11.4000000;
 y1=-17.1000000;y2=-3.0400000;y3=8.5500000;y4=6.2700000;


%% shape functions derivatives in local coordinates

%  row contains derivatives for N1,N2,N3,N4 wrt eta...
for i =1:4
   DN1_eta(i)=   zi(i)/4 - 1/4
   DN2_eta(i)=  -zi(i)/4 - 1/4
   DN3_eta(i)=   zi(i)/4 + 1/4
   DN4_eta(i)=    1/4 - zi(i)/4


% row contains derivatives for N1,N2,N3,N4 wrt zi...

   DN1_zi(i)=((eta(i)/4) - 1/4)
   DN2_zi(i)=( 1/4 - (eta(i)/4)) 
   DN3_zi(i)= eta(i)/4 + 1/4
   DN4_zi(i)= - eta(i)/4 - 1/4
end


K=zeros(8,8)

for i=1:nodes
    
a(i)=0.25*(y1*(zi(i)-1)+y2*(-1-zi(i))+y3*(1+zi(i))+y4*(1-zi(i)))  
b(i)=0.25*(y1*(eta(i)-1)+y2*(1-eta(i))+y3*(1+eta(i))+y4*(-1-eta(i)))    
c(i)=0.25*(x1*(eta(i)-1)+x2*(1-eta(i))+x3*(1+eta(i))+x4*(-1-eta(i))) 
d(i)=0.25*(x1*(zi(i)-1)+x2*(-1-zi(i))+x3*(1+zi(i))+x4*(1-zi(i))) 
    
CC=[0 1-eta(i)  eta(i)-zi(i)  zi(i)-1;
   eta(i)-1 0 zi(i)+1 -zi(i)-eta(i);
   zi(i)-eta(i) -zi(i)-1 0 eta(i)+1 ;
     1-zi(i) zi(i)+eta(i) -eta(i)-1 0]
  
 
 
 J=(1/8)*[x1 x2 x3 x4]*CC*[y1 y2 y3 y4]' 
 
 %% B (derivative of the shape function wrt to guassian points
 B1=[ (a(i)*DN1_zi(i))-(b(i)*DN1_eta(i)) 0;
     0 (c(i)*DN1_eta(i))-(d(i)*DN1_zi(i)) ;
     (c(i)*DN1_eta(i))-(d(i)*DN1_zi(i))   (a(i)*DN1_zi(i))-(b(i)*DN1_eta(i)) ]

 B2=[ (a(i)*DN2_zi(i))-(b(i)*DN2_eta(i)) 0 ;
      0 (c(i)*DN2_eta(i))-(d(i)*DN2_zi(i));
     (c(i)*DN2_eta(i)-d(i)*DN2_zi(i)) (a(i)*DN2_zi(i)-b(i)*DN2_eta(i))]
 B3=[ (a(i)*DN3_zi(i))-(b(i)*DN3_eta(i)) 0;
        0 (c(i)*DN3_eta(i))-(d(i)*DN3_zi(i));
     (c(i)*DN3_eta(i)-d(i)*DN3_zi(i)) (a(i)*DN3_zi(i)-b(i)*DN3_eta(i))]
 B4=[ (a(i)*DN4_zi(i))-(b(i)*DN4_eta(i)) 0;
        0 (c(i)*DN4_eta(i))-(d(i)*DN4_zi(i));
     (c(i)*DN4_eta(i)-d(i)*DN4_zi(i)) (a(i)*DN4_zi(i)-b(i)*DN4_eta(i))]
 B=(1/J)*[B1 B2 B3 B4]

PS=(E/(1-v.^2))*[1 v 0;v 1 0;0 0 0.5*(1-v)] 

K=K+(W1.*W2.*t.*B'*PS*(B).*J)
end
surf(K)

 K_check=[K(1,1),K(2,6),K(4,4),K(8,8)]
 
 mesh(K)
 