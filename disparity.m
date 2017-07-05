function [] = disparity( l_img, r_img, scalefactor, output )

%University of California, Santa Barbara
%Winter 2015
%CS/ECE 181b - Intro to Computer Vision
%Raphael Ruschel dos Santos

occ = 50; %Penalty for the occlusion

left_img = imread(l_img);
right_img = imread(r_img);

left_img = rgb2gray(left_img);
right_img = rgb2gray(right_img);

[M,N,~] = size(left_img); %Both images are the same size

dis_map = zeros(M,N); %Create an empty matrix to put the final result

D0 = zeros(N,N);
A = zeros(N-1,N-1); %I made matrix A with size N-1 because otherwise it 
%would have a row and a column full of zeros and the program would get 
%into an infinite loop during the switch case

D0(1:N,1) = (1:N)*occ; %Populate the original D matrix with the penalties
D0(1,1:N) = (1:N)*occ;

mins = zeros(1,4);

for k=1:M
    D = D0;
    l_img = left_img(k,:); %Extract just one row of each image
    r_img = right_img(k,:);
    for i=2:N
        for j=2:N
            %Calculate the cost to match the pixels
            mins(1) = D(i-1,j-1) + (double(l_img(i)) - double(r_img(j)))^2; 
            mins(2) = D(i-1,j) + occ;
            mins(3) = D(i,j-1) + occ;
            mins(4) = D(i-1,j-1) + 2*occ;
            D(i,j) = min(mins);
            %Set the flags according to the case with the lowest cost
            if(mins(1) == D(i,j)) 
                A(i-1,j-1) = 1; 
            %end;
            elseif(mins(2) == D(i,j)) 
                A(i-1,j-1) = 2; 
            %end;
            elseif(mins(3) == D(i,j)) 
                A(i-1,j-1) = 3; 
            %end;
            %if(mins(4) == D(i,j)) 
            else
                A(i-1,j-1) = 4; 
            end;
        end
    end
    
    
    p = N-1;
    q = N-1;
    
    while(p~=0 && q~=0)
        switch(A(p,q))
            case 1 %p matchs q
                dis_map(k,p) = p - q;
                p = p - 1;
                q = q - 1;
            case 2 %p is unmatched
                p = p - 1;
            case 3 %q is unmatched
                q = q - 1;
            case 4 %both are unmatched
                p = p - 1;
                q = q - 1;
        end
    end
    disp(['Row: ' num2str(k)])
end
dis_map = medfilt2(dis_map);
dis_map = scalefactor*(dis_map/max(dis_map(:)));
imwrite(dis_map,output);

end

