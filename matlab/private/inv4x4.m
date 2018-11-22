function b = inv4x4(a)

b = a;

% accumulating the 4 consecutive tmp-variables yields the 4x4 determinant
tmp = det3x3(a([2 3 4],[2 3 4],:,:));
D   = tmp.*a(1,1,:,:); 
b(1,1,:,:) =  tmp;
tmp = -det3x3(a([2 3 4],[1 3 4],:,:));
D   = D + a(1,2,:,:).*tmp;
b(2,1,:,:) =  tmp;
tmp = det3x3(a([2 3 4],[1 2 4],:,:));
D   = D + a(1,3,:,:).*tmp;
b(3,1,:,:) =  tmp;
tmp = -det3x3(a([2 3 4],[1 2 3],:,:));
D   = D + a(1,4,:,:).*tmp;
b(4,1,:,:) =  tmp;

b(1,2,:,:) = -det3x3(a([1 3 4],[2 3 4],:,:));
b(2,2,:,:) =  det3x3(a([1 3 4],[1 3 4],:,:));
b(3,2,:,:) = -det3x3(a([1 3 4],[1 2 4],:,:));
b(4,2,:,:) =  det3x3(a([1 3 4],[1 2 3],:,:));
b(1,3,:,:) =  det3x3(a([1 2 4],[2 3 4],:,:));
b(2,3,:,:) = -det3x3(a([1 2 4],[1 3 4],:,:));
b(3,3,:,:) =  det3x3(a([1 2 4],[1 2 4],:,:));
b(4,3,:,:) = -det3x3(a([1 2 4],[1 2 3],:,:));
b(1,4,:,:) = -det3x3(a([1 2 3],[2 3 4],:,:));
b(2,4,:,:) =  det3x3(a([1 2 3],[1 3 4],:,:));
b(3,4,:,:) = -det3x3(a([1 2 3],[1 2 4],:,:));
b(4,4,:,:) =  det3x3(a([1 2 3],[1 2 3],:,:));
b          = bsxfun(@rdivide,b,D);
%out = (1:4)';
%in  = [2 3 4;1 3 4;1 2 4;1 2 3];
%a_adj = a;
%s     = [1 -1 1 -1;-1 1 -1 1;1 -1 1 -1;-1 1 -1 1];
%for k = 1:4
%  for m = 1:4
%    a_adj(out(m),out(k),:,:) = s(k,m).*det3x3(a(in(k,:),in(m,:),:,:));
%  end
%end
%b = bsxfun(@rdivide,a_adj,det4x4(a));
%b = a_adj./det4x4(a);

%function b = det3x3(a)
%
%b = a(1,1,:,:).*det2x2(a(2:3,2:3,:,:)) - a(1,2,:,:).*det2x2(a(2:3,[1 3],:,:)) + a(1,3,:,:).*det2x2(a(2:3,1:2,:,:));

%function b = det4x4(a)
%
%b = a(1,1,:,:).*det3x3(a(2:4,2:4,:,:)) - a(1,2,:,:).*det3x3(a(2:4,[1 3 4],:,:)) + a(1,3,:,:).*det3x3(a(2:4,[1 2 4],:,:)) - a(1,4,:,:).*det3x3(a(2:4,1:3,:,:));