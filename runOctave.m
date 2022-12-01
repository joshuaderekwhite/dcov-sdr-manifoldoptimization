X = importdata('X.csv');
Y = importdata('Y.csv');
labels = ostrsplit(cell2mat(X.textdata),",");

[X_sdr, beta, labels_sdr] = ssdr(X.data,Y.data,labels);
