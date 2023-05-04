
close all;
data_test_op = load('test_out.txt');

[num,~] = size(data_test_op);

pCoM = data_test_op(:,1:3);%1
vCoM = data_test_op(:,4:6);%4
peW = data_test_op(:,7:12);%7
legSptInd = data_test_op(:,13:14);%13

L = data_test_op(:,15:17);%15
pf2cen = data_test_op(:,18:20);%18
vw = data_test_op(:,21:23);%21
Lm = data_test_op(:,24:26);%24

psw_Tkf_t = data_test_op(:,27:29);%27
L_Tkf_t = data_test_op(:,30:32);%30
L_des_Tk1f = data_test_op(:,33:35);%33

L_Tkf_t_1 = zeros(num,2);
L_Tkf_t_2 = zeros(num,2);

L_Tkf_t_1(:,1) = data_test_op(:,36);%36
L_Tkf_t_2(:,1) = data_test_op(:,37);%37
L_Tkf_t_1(:,2) = data_test_op(:,38);%38
L_Tkf_t_2(:,2) = data_test_op(:,39);%39

t = data_test_op(:,40);%40
swing = data_test_op(:,41);%41

pcen_w = data_test_op(:,42:44);%42-44

ww = data_test_op(:,45:47);%45-47

Lcom=data_test_op(:,48:50);

figure;%1
for i = 1:3
    subplot(4,3,i);
    plot(pCoM(:,i));
    legend("pCoM");

    subplot(4,3,i+3);
    plot(vCoM(:,i));
    legend("vCoM");

    subplot(4,3,i+6);
    plot(peW(:,i));
    legend("peWr");    

    subplot(4,3,i+9);
    plot(peW(:,i+3));
    legend("peWl");  
end

figure;%2
for i = 1:3
    subplot(4,3,i);
    plot(Lm(:,i),'r');
    hold on;
    plot(L(:,i),'b');
    legend("Lm", "L");

    subplot(4,3,i+3);
    plot(vCoM(:,i));
    legend("vCoM");

    subplot(4,3,i+6);
    plot(ww(:,i));
    legend("ww");

    subplot(4,3,i+9);
    plot(pf2cen(:,i));
    legend("pf2cen");
end

figure;%3
for i = 1:3
    subplot(4,3,i);
    plot(peW(:,i));
    legend("peWr");

    subplot(4,3,i+3);
    plot(peW(:,i+3));
    legend("peWl");

    subplot(4,3,i+6);
    plot(pcen_w(:,i));
    legend("p_cen_w");

    subplot(4,3,i+9);
    plot(pf2cen(:,i));
    legend("pf2cen");

end

figure;%4
for i = 1:3
    subplot(3,3,i);
    plot(psw_Tkf_t(:,i),'b');
    hold on;
%     plot(data_foot_out(:,14),'r');
    legend("psw_Tkf_t");

    subplot(3,3,i+3);
    plot(L_Tkf_t(:,i),'b');
    hold on;
%     plot(data_foot_out(:,14),'r');
    legend("L_Tkf_t");

    subplot(3,3,i+6);
    plot(L_des_Tk1f(:,i),'b');
    hold on;
%     plot(data_foot_out(:,14),'r');
    legend("L_des_Tk1f");
end

figure;%5
subplot(3,2,1);
plot(L_Tkf_t_1(:,1),'b');
hold on;
%     plot(data_foot_out(:,14),'r');
legend("L_Tkf_t_1_x");
subplot(3,2,2);
plot(L_Tkf_t_2(:,1),'b');
hold on;
%     plot(data_foot_out(:,14),'r');
legend("L_Tkf_t_2_x");
subplot(3,2,3);
plot(L_Tkf_t_1(:,2),'b');
hold on;
%     plot(data_foot_out(:,14),'r');
legend("L_Tkf_t_1_y");
subplot(3,2,4);
plot(L_Tkf_t_2(:,2),'b');
hold on;
%     plot(data_foot_out(:,14),'r');
legend("L_Tkf_t_2_y");

subplot(3,2,5);
plot(t,'b');
legend("t");
subplot(3,2,6);
plot(swing,'b');
legend("swing");


figure();
subplot(2,1,1)
plot(Lcom(:,1));
hold on;
plot(L(:,1));
legend('Lcom-x','L-x');
subplot(2,1,2)
plot(Lcom(:,2));
hold on;
plot(L(:,2));
legend('Lcom-y','L-y');







