% Generation of members in the symplectic group Sp(2n)
%
% Author: Michel Barbeau, Carleton University
% Version: November 15, 2018
%
clear;
% usage examples:
%
% disp(mysymplectic(1,2));
% disp(corder(1));disp(corder(2));disp(corder(5));
% disp(sorder(1));disp(sorder(2));disp(sorder(2)/sorder(1));disp(sorder(3));
% disp(sorder(3)/sorder(2));

% [h1, h2] = findxvection([1; 0], [1; 0])
% [h1, h2] = findxvection([1; 0], [1; 1])
% [h1, h2] = findxvection([1; 0; 1; 0; 0; 1; 0; 1], [1; 0; 0; 0; 0; 0; 0; 0])
% findxvectionQA(3);
% v = de2bi(3,2*2,'left-msb')';
% w = de2bi(12,2*2,'left-msb')';
% [h1,h2] = findxvection(v, w);
% t1=xvection(h1,v);
% t2=xvection(h2,w);
% t1~=t2
% [v, w, h1, h2, t1, t2]à
% [v, w, t, h0, h1, h2] = indextosp(48,2); disp([v, w, t, h0, h1, h2]);
% indextospQA(2);
% sigma = mysymplectic(48,2); disp(sigma);
mysymplecticQA(2); % greater than 2 takes a long time

function [m] = sorder(n)
% Returns the order of symplectic group Sp(2n)
m = 1;
for j=1:n
    m = m * (power(4,j) - 1);
end
m = power(2,power(n,2)) * m;
end

function [m] = porder(n)
% Returns the order of group "S_n"
if n==1
    m = sorder(1);
else
    m = sorder(n)/sorder(n-1);
end
end

function [v,w,s,h0,h1,h2] = indextosp(i,n)
% Given an index "i" and dimension "n", returns the
% i-th symplectic pair in group "S_n"
% index validation
if i<0 || i>=porder(n)
	error('symplectic pair index out of range n=%d i=%d', n, i)
end
% decompose order into two factors
d = power(2,2*n-1);
r = floor(i/d) + 1;
% "v" becomes the binary expansion of "r" 
v = de2bi(r,2*n,'left-msb')';
% find transvections that yield "v" from "I_{*,1}"
[h1,h2] = findxvection([ 1 zeros(1,length(v)-1)]',v);
% decompose "i" into 2nd part
s = mod(i,d);
% binary expansion of "t" 
b = de2bi(s,2*n-1,'left-msb');
% apply xvections to I_{*,1}+b(2)I_{*,3}+...+b(2*n-1)I_{*,2n}
I = eye(2*n); % identity matrix of dimension 2n 
h0 = xvection(h1,xvection(h2, sum(I(:,find([1 0 b(2:end)])),2) ));  
% construct transvection "t"
s = ~b(1)*v;
% generate "w"
w =  xvection(s,xvection(h0,xvection(h1,...
    xvection(h2,[ 0 1 zeros(1,length(v)-2)]'))));
end

function [L] = Lambda(n) 
% returns the Lambda(n) matrix
    L = kron(eye(n),[0 1; 1 0]);
end

function [p] = sip(v,w) 
% symplectic inner product
   p = mod(v'*Lambda(length(v)/2)*w,2);
end

function [w] = xvection(h,v)
% Applies the transvection "Z_h" to column vector "v"
   w = mod(v + sip(v,h)*h,2);
end

function [h1,h2] = findxvection(v,w)
% Gigen two colum vectors "v" and "w",
% find h1 and h2 such that w = Z_h1 Z_h2 w
%
% initialize transvection vectors
h1 = zeros(length(v),1);
h2 = zeros(length(v),1);
if isequal(v,w) % vectors equal?
    return;
end
if sip(v,w) % symplectic inner product is one? 
    h1 = xor(v,w);
    return;
end
z = zeros(length(v),1);
% find v~=00 and w~=00
for j=2:2:length(v)
    if (v(j-1)||v(j)) && (w(j-1)||w(j))
        % pair is found!
        z(j-1) = xor(v(j-1),w(j-1));
        z(j) = xor(v(j),w(j));
        if (z(j-1)+z(j))==0
            z(j)=1;
            if v(j-1)~=v(j)
                z(j-1)=1;
            end
         end
         h1 = xor(v,z);
         h2 = xor(w,z);
         return;
    end
end
% find v~=00 and w=00
for j=2:2:length(v)
    if (v(j-1)||v(j)) && ~(w(j-1)||w(j))
        % pair is found
        if v(j-1)==v(j)
            z(j)=1;
        else
            z(j)=v(j-1);
            z(j-1)=v(j);
        end
        break;
    end
end
% find v==00 and w~=00
for j=2:2:length(v)
    if ~(v(j-1)||v(j)) && (w(j-1)||w(j))
        % pair is found
        if w(j-1)==w(j)
            z(j)=1;
        else
            z(j)=w(j-1);
            z(j-1)=w(j);
        end
        break;
    end
end
h1 = xor(v,z);
h2 = xor(w,z);
end

function [m3] = directsum(m1,m2)
% Input:
%    m1, m2 = two square matrices
% Output: m3 equal to
%    ( m1 0  )
%    ( 0  m2 )
% Create a square matrix of zeros, dimension is
% the sum of dimensions of m1 and m2
m3 = vertcat( ...
   horzcat(m1,zeros(length(m1),length(m2))),...
   horzcat(zeros(length(m2),length(m1)), m2) );
end

function [sigma] = mysymplectic(i,n)
% Given an index "i" and dimension "n", returns the
% i-th symplectic in group Sp(2n).
% 
% validation of inputs
if mod(i,1) || i<0
    error('1st argument must be non negative integer i=%d', i)
end
if mod(i,1) && n<=0
    error('2nd argument must be a positive integer n=%d', n)
end
if i<0 || i>=sorder(n) % symplectics indexed in 0...|Sp(2n)|-1
    error('symplectic index out of range n=%d i=%d', n, i)
end
% order of symplectic pair group "S_n"
d = porder(n);
% find symplectic pair and transvections for index "i"
[v, w, t, h0, h1, h2] = indextosp(mod(i,d),n);
if n==1 % done!
    sigma = [ v w ];
    return;
else
    sigma = directsum(eye(2),mysymplectic(floor(i/d),n-1));
    for i=1:2*n
        % apply transvections to column "i"
        sigma(:,i) = xvection(t,xvection(h0,xvection(h1,...
          xvection(h2,sigma(:,i)))));
    end
end
end

% Qualiy Assurance Code

function [m] = corder(n)
% Returns the order of the Cliffor group C(n)
m = 1;
for j=1:n
    m = m * (power(4,j) - 1);
end
m = power(2,power(n,2)+2*n) * m;
end

function findxvectionQA(n)
% Verifies function findxvection for all non null binary column vectors of 
% length "n"
% max value on a vector of length 2n
max = power(2,2*n)-1;
for i=1:max
    for j=1:max
        v = de2bi(i,2*n,'left-msb')';
        w = de2bi(j,2*n,'left-msb')';
        [h1,h2] = findxvection(v, w);
        if ~isequal(w,xvection(h1,xvection(h2,v)))
            disp(v); disp(w); disp(h1); disp(h2);
            error('findxvection errror i=%d j=%d', i, j);      
        end
    end
end
end

function indextospQA(n)
% Verifies function indextosp
for i=0:porder(n)-1
    [v, w, t, h0, h1, h2] = indextosp(i,n);
    % disp([v, w, t, h0, h1, h2]);
    if ~sip(v,w)
        disp([v, w, t, h0, h1, h2]);
        error('indextosp errror n=%d i=%d', n, i);      
    end
end
end

function mysymplecticQA(n)
% Verifies function mysymplectic
% generate the Lambda(n) matrix
L = Lambda(n);
% order of Sp(2n)
d =  sorder(n);
% exhaustive verification
fprintf('Generating %d symplectics\n', d);
for i=0:d-1
    sigma(:,:,i+1) = mysymplectic(i,n);
    if ~isequal(mod(sigma(:,:,i+1)*L*sigma(:,:,i+1)',2),L)
        disp(sigma(:,:,i+1));
        error('symplectic test failed n=%d i=%d', n, i);      
    end
end
fprintf('%d symplectics generated\n', d);
% verification of uniqueness
for i=0:d-1
    if i<d-1
        % compare with every othe matrics from that indx
        for j=i+1:d-1
            if isequal(sigma(:,:,i+1),sigma(:,:,j+1))
               disp(sigma(:,:,i+1),sigma(:,:,j+1));
               error('repetition n=%d i=%d j=%d',n, i,j);
            end
        end
    end
end
fprintf('uniqueness test passed\n')
end
