function [ret_val] = wedge_block(w,c)
a=ones(w,w);
b=ones(w,w);
%% case1
if c==1
    for i=1:w
        for j=1:w
            if i<=ceil(w/2)
                if j<=ceil(w/2)
                    if(j>=i)
                        a(i,j)=0.5;
                    end
                else
                    if(((j-ceil(w/2))+i)<=(ceil(w/2)+1))
                        a(i,j)=0.5;
                    end
                end
            end
        end
    end
    
    for i=1:w
        for j=1:w
            if i>ceil(w/2)
                if j>=ceil(w/2)
                    if(j<=i)
                        b(i,j)=0.5;
                    end
                else
                    if(((j-ceil(w/2))+i)>=(ceil(w/2)+1))
                        b(i,j)=0.5;
                    end
                end
            else
            end
        end
    end
else
    % imagesc(a)
    % imagesc(a);
    
    for i=1:w
        for j=1:w
            if j<=ceil(w/2)
                if i<=ceil(w/2)
                    if(j<=i)
                        a(i,j)=0.5;
                    end
                else
                    if(((i-ceil(w/2))+j)<=(ceil(w/2)+1))
                        a(i,j)=0.5;
                    end
                end
            else
            end
        end
    end
    % imagesc(a);
    for i=1:w
        for j=1:w
            if j>=ceil(w/2)
                if i<=ceil(w/2)
                    if((j-ceil(w/2)+i)>=(ceil(w/2)+1))
                        b(i,j)=0.5;
                    end
                else
                    if(j>=i)
                        b(i,j)=0.5;
                    end
                end
            end
        end
    end
    % imagesc(a);
    
end
[px, py] = meshgrid(linspace(-1, 1, w));
x=px(floor(w/2),:);
y=py(floor(w/2),:);
pxy=sqrt(x.^2 + y.^2);
smooth_edge=ones(w,w);
k=floor(w/2);
for i=1:w
    if i<=floor(w/2)
        smooth_edge(i,:)=(1 + erf(10 * (pxy - (w-(i-1)*2) / w)));
    else
        smooth_edge(i,:)= smooth_edge(k,:);
        k=k-1;
    end
end
if c==1
    ret_val = a.*b.*smooth_edge-0.1;
else
    ret_val = a.*b.*smooth_edge'-0.1;
end

end


