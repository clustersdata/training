# training
training
clear all;
close all;
clc;

n=5;
randn('seed',1);
mu1=[0 0];
S1=[0.5 0;
    0 0.5];
P1=mvnrnd(mu1,S1,n);

mu2=[0 6];
S2=[0.5 0;
    0 0.5];
P2=mvnrnd(mu2,S2,n);

mu3=[6 3];
S3=[0.5 0;
    0 0.5];
P3=mvnrnd(mu3,S3,n);


P=[P1;P2;P3];
meanP=mean(P);

P=[P(:,1)-meanP(1) P(:,2)-meanP(2)];

sigma = 5;

X=P(:,1);
Y=P(:,2);
B=rand(3*n,1);

w1 = rand(3*n,1);
w2 = rand(3*n,1);
w3 = rand(3*n,1);

w4 = rand(3*n,1);
w5 = rand(3*n,1);
w6 = rand(3*n,1);


for i=1:3*n
    i
    while 1
        
        y1 = X(i)*w1(i) + Y(i)*w4(i) + B(i);       
        y2 = X(i)*w2(i) + Y(i)*w5(i) + B(i);        
        y3 = X(i)*w3(i) + Y(i)*w6(i) + B(i);     
        
        h1 = 1/(1+exp(-y1));
        h2 = 1/(1+exp(-y2));       
        h3 = 1/(1+exp(-y3));      
        
        e1  = 1/2*(1 - h1)^2;
        e2  = 1/2*(1 - h2)^2;       
        e3  = 1/2*(1 - h3)^2;
 
        if i<=n && e1<=0.0000001
            break;
        elseif i>n && i<=2*n && e2<0.0000001
            break;
        elseif i>2*n && e3<0.0000001
            break;
        end
        
        
        if i<=n
            w1(i) = w1(i)-sigma*(h1-1)*h1*(1-h1)*X(i);
            w2(i) = w2(i)-sigma*(h2-0)*h2*(1-h2)*X(i);
            w3(i) = w3(i)-sigma*(h3-0)*h3*(1-h3)*X(i);    
            
            w4(i) = w4(i)-sigma*(h1-1)*h1*(1-h1)*Y(i);
            w5(i) = w5(i)-sigma*(h2-0)*h2*(1-h2)*Y(i);
            w6(i) = w6(i)-sigma*(h3-0)*h3*(1-h3)*Y(i);                   
            
            B(i) =B(i)- sigma*((h1-1)*h1*(1-h1)+(h2-0)*h2*(1-h2)+(h3-0)*h3*(1-h3));
        elseif i>n && i<=2*n
            w1(i) = w1(i)-sigma*(h1-0)*h1*(1-h1)*X(i);
            w2(i) = w2(i)-sigma*(h2-1)*h2*(1-h2)*X(i);
            w3(i) = w3(i)-sigma*(h3-0)*h3*(1-h3)*X(i);    
            
            w4(i) = w4(i)-sigma*(h1-0)*h1*(1-h1)*Y(i);
            w5(i) = w5(i)-sigma*(h2-1)*h2*(1-h2)*Y(i);
            w6(i) = w6(i)-sigma*(h3-0)*h3*(1-h3)*Y(i);                   
            
            B(i) =B(i)- sigma*((h1-0)*h1*(1-h1)+(h2-1)*h2*(1-h2)+(h3-0)*h3*(1-h3));         
        else
            w1(i) = w1(i)-sigma*(h1-0)*h1*(1-h1)*X(i);
            w2(i) = w2(i)-sigma*(h2-0)*h2*(1-h2)*X(i);
            w3(i) = w3(i)-sigma*(h3-1)*h3*(1-h3)*X(i);    
            
            w4(i) = w4(i)-sigma*(h1-0)*h1*(1-h1)*Y(i);
            w5(i) = w5(i)-sigma*(h2-0)*h2*(1-h2)*Y(i);
            w6(i) = w6(i)-sigma*(h3-1)*h3*(1-h3)*Y(i);                   
            
            B(i) =B(i)- sigma*((h1-0)*h1*(1-h1)+(h2-0)*h2*(1-h2)+(h3-1)*h3*(1-h3));                   
        end
         

    end
end

plot(P(:,1),P(:,2),'o');
hold on;

flag = 0;
M=[];
for x=-8:0.3:8
    for y=-8:0.3:8

        H=[]; 
        for i=1:3*n
            y1 = x*w1(i)+y*w4(i) +B(i);
            y2 = x*w2(i)+y*w5(i) +B(i);
            y3 = x*w3(i)+y*w6(i) +B(i);
            h1=1/(1+exp(-y1));
            h2=1/(1+exp(-y2));
            h3=1/(1+exp(-y3));
            
            H=[H;h1 h2 h3];
        end
  %      H1 = mean(H(1:n,1));
  %      H2 = mean(H(n:2*n,2));
  %      H3 = mean(H(2*n:3*n,3));
        
        meanH = mean(H);
        H1 = meanH(1);
        H2 = meanH(2);
        H3= meanH(3);
        if H1>H2 && H1>H3
            plot(x,y,'g.')
        elseif H2 > H1 && H2 > H3
            plot(x,y,'r.')
        elseif H3 > H1 && H3 > H2
            plot(x,y,'b.')
        end
        
    end
end
