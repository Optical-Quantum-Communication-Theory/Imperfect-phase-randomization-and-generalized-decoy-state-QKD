clear;
results = load('results.mat');

keyRate = results.keyRate;
distance = results.distance;
probPhaseRand = results.a;

semilogy(distance,keyRate(1,:),'-o',"DisplayName", sprintf("q = %f",probPhaseRand(1)));
hold on
semilogy(distance,keyRate(2,:),'-+',"DisplayName", sprintf("q = %f",probPhaseRand(2)));
semilogy(distance,keyRate(3,:),'-*',"DisplayName", sprintf("q = %f",probPhaseRand(3)));
hold off
xlabel('Distance (km)')
ylabel('Key Rate')
legend()