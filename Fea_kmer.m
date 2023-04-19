function kmer_Matrix=Fea_kmer(file_name,k)
[hp,Seq]=fastaread(file_name);
Np=length(hp); %number of samples
%% kmer
AA='ACGT';
m=length(AA);%
r=0;
H=cell(1,m^k);
Wp=zeros(Np,m^k);

switch k
    case {1}
        for i=1:m
            x1=AA(i);
            r=r+1;
            H{1,r}=x1; 
        end
        for i=1:Np
            for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a=s(j);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
            end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
        end  
        
    case {2}
        for i=1:m
           x1=AA(i);
            for j=1:m
                x2=AA(j);
                x=strcat(x1,x2);
                r=r+1;
                H{1,r}=x;           
            end  
        end
        for i=1:Np
            for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a1=s(j);
            a2=s(j+1);
            a=strcat(a1,a2);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
            end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
        end  
    case {3}
     for i=1:m
           x1=AA(i);
        for j=1:m
             x2=AA(j);
            for t=1:m
               x3=AA(t);
               x=strcat(x1,x2,x3);
               r=r+1;
               H{1,r}=x;           
            end
        end  
     end
     for i=1:Np
        for j=1:(length(Seq{1,i})-k+1)
        s=Seq{1,i};
        a1=s(j);
        a2=s(j+1);
        a3=s(j+2);
        a=strcat(a1,a2,a3);
        g=strmatch(a,H,'exact');
        Wp(i,g)=Wp(i,g)+1;
        end
        Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
     end  
     
     case {4}
     for i=1:m
           x1=AA(i);
        for j=1:m
             x2=AA(j);
            for t=1:m
               x3=AA(t);
               for e=1:m
                   x4=AA(e);
                   x=strcat(x1,x2,x3,x4);
                   r=r+1;
                   H{1,r}=x; 
               end
            end
        end  
     end
     for i=1:Np
         for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a1=s(j);
            a2=s(j+1);
            a3=s(j+2);
            a4=s(j+3);
            a=strcat(a1,a2,a3,a4);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
         end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
     end  
    case {5}
     for i=1:m
           x1=AA(i);
        for j=1:m
             x2=AA(j);
            for t=1:m
               x3=AA(t);
               for e=1:m
                   x4=AA(e);
                   for e1=1:m
                       x5=AA(e1);
                   x=strcat(x1,x2,x3,x4,x5);
                   r=r+1;
                   H{1,r}=x; 
                   end
               end
            end
        end  
     end
     for i=1:Np
         for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a1=s(j);
            a2=s(j+1);
            a3=s(j+2);
            a4=s(j+3);
            a5=s(j+4);
            a=strcat(a1,a2,a3,a4,a5);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
         end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
     end  
   case {6}
     for i=1:m
           x1=AA(i);
        for j=1:m
             x2=AA(j);
            for t=1:m
               x3=AA(t);
               for e=1:m
                   x4=AA(e);
                   for e1=1:m
                       x5=AA(e1);
                       for e2=1:m
                           x6=AA(e2);
                   x=strcat(x1,x2,x3,x4,x5,x6);
                   r=r+1;
                   H{1,r}=x; 
                       end
                   end
               end
            end
        end  
     end
     for i=1:Np
         for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a1=s(j);
            a2=s(j+1);
            a3=s(j+2);
            a4=s(j+3);
            a5=s(j+4);
            a6=s(j+5);
            a=strcat(a1,a2,a3,a4,a5,a6);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
         end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
     end  
  case {7}
     for i=1:m
           x1=AA(i);
        for j=1:m
             x2=AA(j);
            for t=1:m
               x3=AA(t);
               for e=1:m
                   x4=AA(e);
                   for e1=1:m
                       x5=AA(e1);
                       for e2=1:m
                           x6=AA(e2);
                           for e3=1:m
                               x7=AA(e3);
                   x=strcat(x1,x2,x3,x4,x5,x6,x7);
                   r=r+1;
                   H{1,r}=x; 
                           end
                       end
                   end
               end
            end
        end  
     end
     for i=1:Np
         for j=1:(length(Seq{1,i})-k+1)
            s=Seq{1,i};
            a1=s(j);
            a2=s(j+1);
            a3=s(j+2);
            a4=s(j+3);
            a5=s(j+4);
            a6=s(j+5);
            a7=s(j+6);
            a=strcat(a1,a2,a3,a4,a5,a6,a7);
            g=strmatch(a,H,'exact');
            Wp(i,g)=Wp(i,g)+1;
         end
            Wp(i,:)=Wp(i,:)/(length(Seq{1,i})-k+1);
     end  
     otherwise
        warning('Unexpected k value')
end
%%
kmer_Matrix=Wp;
 
 
