load('Psi_eyx_xn.mat');
aa=Psi_eyx_xn;
a=cat(3,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,...
    a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,...
    a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,...
    a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,...
    a40,a41,a42);

b=(aa-a);
subplot(1,3,1);
surf(a(:,:,13));
subplot(1,3,2);
surf(aa(:,:,13));
subplot(1,3,3);
plot(reshape(b,[1,prod(size(aa))]));
