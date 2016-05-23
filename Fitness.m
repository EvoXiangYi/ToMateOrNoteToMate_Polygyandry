function [fitnessVector]=Fitness(habitatVector, envTraitVector, habitatCenters, fitnessIOS)
	popSize=length(habitatVector);
	numHabitats=length(habitatCenters);
	fitnessVector=zeros(1,popSize);
	for i=1:numHabitats
		idx = find(habitatVector==i);
		distances = abs(envTraitVector(idx) - habitatCenters(i));
		fitness = exp(- fitnessIOS * distances .^ 2);
		fitnessVector(idx) = fitness;
	end
end


%%%%%%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%
%habitatVector=[1,2,3,4,1,2,3,4,1,2,3,4,0,0];
%envTraitVector=[0.05,0.15,0.25,0.35,0.1,0.2,0.3,0.4,0.15,0.25,0.35,0.45,0.5,1];
%habitatCenters=[0.1,0.2,0.3,0.4];
%fitnessIOS=3;
%fitnessVector=Fitness(habitatVector, envTraitVector, habitatCenters, fitnessIOS)