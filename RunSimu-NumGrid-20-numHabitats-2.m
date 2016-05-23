numRepeats=10;
for i=1:numRepeats
	x=[3000,20,1000,2,0.9,0.9,3,0,1,0.01,0.01,0.01];
	x=mat2cell(x,1,ones(1,length(x)));
	x{length(x)+1}=i;
	DispersalSimulation_mateOrNot_Polygyandry(x{:});
end