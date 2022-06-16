clc;
% 目标域数据划分为训练数据和测试数据，在求解Pim、P、Zim等变量过程中只使用目标域训练数据和源域所有数据
% 获取目标域训练数据
% class 1
load Target_ASD_train.txt; % shape：[样本数，样本特征维度]
Tm = Target_ASD_train';
% class 2
load Target_NC_train.txt;
Tn = Target_NC_train';

% 获取目标域测试数据
load Target_test.txt;
T_test = Target_test';

% 获取源域1数据
% class 1
load Source1_ASD.txt;
S1m = Source1_ASD;
% class 2
load Source1_NC.txt;
S1n = Source1_NC.txt;

% 获取源域2数据
% class 1
load Source2_ASD.txt;
S2m = Source2_ASD;
% class 2
load Source2_NC.txt;
S2n = Source2_NC.txt;

% 获取源域3数据
% class 1
load Source3_ASD.txt;
S3m = Source3_ASD;
% class 2
load Source3_NC.txt;
S3n = Source3_NC.txt;

% 获取源域4数据
% class 1
load Source4_ASD.txt;
S4m = Source4_ASD;
% class 2
load Source4_NC.txt;
S4n = Source4_NC.txt;

% 获取源域5数据
% class 1
load Source5_ASD.txt;
S5m = Source5_ASD;
% class 2
load Source5_NC.txt;
S5n = Source5_NC.txt;

% 获取源域6数据
% class 1
load Source6_ASD.txt;
S6m = Source6_ASD;
% class 2
load Source6_NC.txt;
S6n = Source6_NC.txt;

Ssm.S1m = S1m;Ssm.S2m = S2m;Ssm.S3m = S3m;Ssm.S4m = S4m;Ssm.S5m = S5m;Ssm.S6m = S6m;
Ssn.S1n = S1n;Ssn.S2n = S2n;Ssn.S3n = S3n;Ssn.S4n = S4n;Ssn.S5n = S5n;Ssn.S6n = S6n;

% 给参数进行赋值
a = 0.0021003;
b = 1.6111e-06;
opts.tol = 1e-7;
opts.max_iter = 120;
opts.rho = 1.1005;
opts.mu = 1.0058e-07;
opts.max_mu = 1e7;
opts.DEBUG = 1;

cd C:\Users\11944\Desktop\CDLRR\CDLRR % 进入CDLRR文件所在路径
[F_sm,F_sn,Z_sm,Z_sn,P_sm,P_sn,P,Es_sm,Es_sn,obj,err] = CDLRR(Tm,Tn,Ssm,Ssn,a,b,opts);

% New Target domain train data (训练数据)
% class 1 
N_Target_ASD_train = P*Tm;
N_Target_ASD_train = N_Target_ASD_train';
fid1 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Target_ASD_train.txt', 'wt');   % 新数据写入路径
[m1, n1] = size(N_Target_ASD_train);
for i1=1:1:m1
    for j1=1:1:n1
        fprintf(fid1,'%f ',N_Target_ASD_train(i1,j1));
    end
    fprintf(fid1,'\n');
end
fclose(fid1);
% class 2
N_Target_NC_train = P*Tn;
N_Target_NC_train =N_Target_NC_train';
fid2 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Target_NC_train.txt', 'wt');
[m2, n2] = size(N_Target_NC_train);
for i2=1:1:m2
    for j2=1:1:n2
        fprintf(fid2,'%f ',N_Target_NC_train(i2,j2));
    end
    fprintf(fid2,'\n');
end
fclose(fid2);

% New Source domains
% Source 1
% class 1
N_Source1_ASD = P_sm.P1m*S1m;
N_Source1_ASD = N_Source1_ASD';
fid3 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source1_ASD.txt', 'wt');
[m3,n3] = size(N_Source1_ASD);
for i3=1:1:m3
    for j3=1:1:n3
        fprintf(fid3,'%f ',N_Source1_ASD(i3,j3));
    end
    fprintf(fid3,'\n');
end
fclose(fid3);
% class 2
N_Source1_NC = P_sn.P1n*S1n;
N_Source1_NC = N_Source1_NC';
fid4 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source1_NC.txt', 'wt');
[m4,n4] = size(N_Source1_NC);
for i4=1:1:m4
    for j4=1:1:n4
        fprintf(fid4,'%f ',N_Source1_NC(i4,j4));
    end
    fprintf(fid4,'\n');
end
fclose(fid4);
% Source 2
% class 1
N_Source2_ASD = P_sm.P2m*S2m;
N_Source2_ASD = N_Source2_ASD';
fid5 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source2_ASD.txt', 'wt');
[m5,n5] = size(N_Source2_ASD);
for i5=1:1:m5
    for j5=1:1:n5
        fprintf(fid5,'%f ',N_Source2_ASD(i5,j5));
    end
    fprintf(fid5,'\n');
end
fclose(fid5);
% class 2
N_Source2_NC = P_sn.P2n*S2n;
N_Source2_NC = N_Source2_NC';
fid6 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source2_NC.txt', 'wt');
[m6,n6] = size(N_Source2_NC);
for i6=1:1:m6
    for j6=1:1:n6
        fprintf(fid6,'%f ',N_Source2_NC(i6,j6));
    end
    fprintf(fid6,'\n');
end
fclose(fid6);
% Source 3
% class 1
N_Source3_ASD = P_sm.P3m*S3m;
N_Source3_ASD = N_Source3_ASD';
fid7 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source3_ASD.txt', 'wt');
[m7,n7] = size(N_Source3_ASD);
for i7=1:1:m7
    for j7=1:1:n7
        fprintf(fid7,'%f ',N_Source3_ASD(i7,j7));
    end
    fprintf(fid7,'\n');
end
fclose(fid7);
% class 2
N_Source3_NC = P_sn.P3n*S3n;
N_Source3_NC = N_Source3_NC';
fid8 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source3_NC.txt', 'wt');
[m8,n8] = size(N_Source3_NC);
for i8=1:1:m8
    for j8=1:1:n8
        fprintf(fid8,'%f ',N_Source3_NC(i8,j8));
    end
    fprintf(fid8,'\n');
end
fclose(fid8);
% Source 4
% class 1
N_Source4_ASD = P_sm.P4m*S4m;
N_Source4_ASD = N_Source4_ASD';
fid9 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source4_ASD.txt', 'wt');
[m9,n9] = size(N_Source4_ASD);
for i9=1:1:m9
    for j9=1:1:n9
        fprintf(fid9,'%f ',N_Source4_ASD(i9,j9));
    end
    fprintf(fid9,'\n');
end
fclose(fid9);
% class 2
N_Source4_NC = P_sn.P4n*S4n;
N_Source4_NC = N_Source4_NC';
fid10 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source4_NC.txt', 'wt');
[m10,n10] = size(N_Source4_NC);
for i10=1:1:m10
    for j10=1:1:n10
        fprintf(fid10,'%f ',N_Source4_NC(i10,j10));
    end
    fprintf(fid10,'\n');
end
fclose(fid10);
% Source 5 
% class 1
N_Source5_ASD = P_sm.P5m*S5m;
N_Source5_ASD = N_Source5_ASD';
fid11 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source5_ASD.txt', 'wt');
[m11,n11] = size(N_Source5_ASD);
for i11=1:1:m11
    for j11=1:1:n11
        fprintf(fid11,'%f ',N_Source5_ASD(i11,j11));
    end
    fprintf(fid11,'\n');
end
fclose(fid11);
% class 2
N_Source5_NC = P_sn.P5n*S5n;
N_Source5_NC = N_Source5_NC';
fid12 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source5_NC.txt', 'wt');
[m12,n12] = size(N_Source5_NC);
for i12=1:1:m12
    for j12=1:1:n12
        fprintf(fid12,'%f ',N_Source5_NC(i12,j12));
    end
    fprintf(fid12,'\n');
end
fclose(fid12);
% Source 6
% class 1
N_Source6_ASD = P_sm.P6m*S6m;
N_Source6_ASD = N_Source6_ASD';
fid13 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source6_ASD.txt', 'wt');
[m13,n13] = size(N_Source6_ASD);
for i13=1:1:m13
    for j13=1:1:n13
        fprintf(fid13,'%f ',N_Source6_ASD(i13,j13));
    end
    fprintf(fid13,'\n');
end
fclose(fid13);
% class 2
N_Source6_NC = P_sn.P6n*S6n;
N_Source6_NC = N_Source6_NC';
fid14 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Source6_NC.txt', 'wt');
[m14,n14] = size(N_Source6_NC);
for i14=1:1:m14
    for j14=1:1:n14
        fprintf(fid14,'%f ',N_Source6_NC(i14,j14));
    end
    fprintf(fid14,'\n');
end
fclose(fid14);

% New Target domian test data (测试集数据)
N_Target_test = P*T_test;
N_Target_test = N_Target_test';
fid15 = fopen('C:\Users\11944\Desktop\CDLRR\CDLRR\New_data\N_Target_test.txt', 'wt');
[m15,n15] = size(N_Target_test);
for i15=1:1:m15
    for j15=1:1:n15
        fprintf(fid15,'%f ',N_Target_test(i15,j15));
    end
    fprintf(fid15,'\n');
end
fclose(fid15);