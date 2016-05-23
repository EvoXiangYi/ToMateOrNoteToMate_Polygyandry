function []= DispersalSimulation_mateOrNot_Polygyandry(InTimeSteps, InNumGrid, InPopSize, InNumHabitats, InpSpatial, InpTemporal, InfitnessIOS, InMeanDispersalTrait_init, InSigmaDispersalTrait_init, InSigmaTraitCoordinates, InSigmaDispersal, InDispersalDeathRate, InNumRepeats)
	
	%%%%%%%%%%%%%%%%%%%%% start Parameters Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%
	numRepeats = InNumRepeats;
	timeSteps = InTimeSteps;
	numGrid = InNumGrid;
	popSize = InPopSize;
	numHabitats = InNumHabitats;
	carryingCapacity = 2 * popSize/numHabitats; %carrying capacity of each of the habitats
	pSpatial=InpSpatial; % Spatial autocorrelation parameter
	pTemporal=InpTemporal; % Temporal autocorrelation parameter
	fitnessIOS=InfitnessIOS;
	init_Mean_Dispersal = InMeanDispersalTrait_init;
	init_Prob_dispersal=1/(1+exp(-init_Mean_Dispersal)); %probability of dispersal for both males and females
	init_Sigma_Dispersal = InSigmaDispersalTrait_init;
	sigmaTraitCoordinates=InSigmaTraitCoordinates;
	sigmaDispersal=InSigmaDispersal;
	dispersalDeathRate=InDispersalDeathRate;
	%%%%%%%%%%%%%%%%%%%%%  end Parameters Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start Habitats Initialization  %%%%%%%%%%%%%%%%%%%%
	%To always start from the same environment
	randomStete=csvread("randStateFile.csv");
	rand("state",randomStete);
	inputEnvironment=rand(1,numGrid); 	%In Luke and Hanna's paper: inputEnvironment=unidrnd(2,1,numGrid)-1;
	envrionmentVector=UpdateEnvironment(inputEnvironment,numGrid,pSpatial,pTemporal,500);
	habitatLocations=SetHabitats(numGrid, numHabitats);
	habitatCenters=envrionmentVector(habitatLocations);
	%%%%%%%%%%%%%%%%%%%%%  end Habitats Initialization  %%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start GenomeStructure  %%%%%%%%%%%%%%%%%%%%%%
	envTraitLocus=1:2; %adaptation traits from both mom and dad
	dispLociF=3:4;
	dispLociS=5:6;
	dispLociM=7:8;
	sex=9;
	habitat=10;
	probMating = 11;
	numLoci=probMating;
	%%%%%%%%%%%%%%%%%%%%%  end GenomeStructure  %%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start PopulationStructure  %%%%%%%%%%%%%%%%%%%%
	population=zeros(numLoci, popSize); %Initialize the "diploid genome"
	population(envTraitLocus,:) = rand(envTraitLocus(end), popSize); %Initialize environmatal coordinates
	population(dispLociF(1):dispLociM(end),:) = normrnd(init_Mean_Dispersal,init_Sigma_Dispersal,[dispLociM(end)-envTraitLocus(end),popSize]);
	population(sex,:) = unidrnd(2,1,popSize)-1; %female:0, males:1
	population(habitat,:) = unidrnd(numHabitats,1,popSize);
	population(probMating,:) = 1;
	%%%%%%%%%%%%%%%%%%%%%  end PopulationStructure  %%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start Output Data Structure  %%%%%%%%%%%%%%%%%%%%
	% 1->Total number of individuals
	% 2->Total numbers of Habitats that have both males and sexual females
	% 3->Mean fitness for females
	% 4->Mean fitness for males
	% 5->Mean fitness for asexual females
	% 6->Mean dispersal probability for sexual females
	% 7->Mean dispersal probability for males
	% 8->Mean dispersal probability for asexual females
	% 9->frequency of asexual females in the whole population
	SummaryData=NaN(9,timeSteps);
	%%%%%%%%%%%%%%%%%%%%%  end Output Data Structure  %%%%%%%%%%%%%%%%%%%% 
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%% start Population Update    %%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for t=1:timeSteps
		t
		
		%%%%%%%%%%%%%%%%%%%%%  start Environment Selection  %%%%%%%%%%%%%%%%%%
	
		%%%%%%%%%%%%%%%%%%%%%  end Environment Selection  %%%%%%%%%%%%%%%%%%%%
	
	
		%%%%%%%%%%%%%%%%%%%%%  start Carrying Capacity regulation  %%%%%%%%%%%%%%%%%%
		surviveVectorAfterRegulation = CarryingCapacityRegulation(population(habitat,:), numHabitats, carryingCapacity);
		population(habitat,surviveVectorAfterRegulation==0)=0;
		%%%%%%%%%%%%%%%%%%%%%  end Carrying Capacity regulation  %%%%%%%%%%%%%%%%%%%%
	
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  start population census  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% 1 ->Total number of individuals : totalNumIndividuals
		totalNumIndividuals=sum(population(habitat,:)!=0);
		if totalNumIndividuals < 1
			break; %exit simulation upon extinction.
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start finding asexual females %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% 9 ->frequency of asexual females in the whole population: freqAsexualFemales
		asexual_females = find(and(rand(1,popSize)>population(probMating,:),population(sex,:)==0));
		population(sex, asexual_females) = 2; %Assign asexual females a different sex.
		numAsexualFemales = length(asexual_females);
		freqAsexualFemales = numAsexualFemales/totalNumIndividuals;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end finding asexual females %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ...
		
		%% 2 ->Total numbers of Habitats that have both males and sexual females : totalEffectiveHabitats
		effectiveHabitats=EffectiveHabitats(numHabitats, population(habitat,:),population(sex,:)); 
		totalEffectiveHabitats=length(effectiveHabitats);
		
		%% 3,4,5 ->Mean fitness for females, males and asexual females : meanFemaleFitness, meanMaleFitness, meanAsexualFemaleFitness
		femalesOnHabitatsCell=cell(1,totalEffectiveHabitats);
		femalesFitnessOnHabitatsCell=cell(1,totalEffectiveHabitats);
		malesOnHabitatsCell=cell(1,totalEffectiveHabitats);
		malesFitnessOnHabitatsCell=cell(1,totalEffectiveHabitats);
		
		for i=1:totalEffectiveHabitats
			localHabitat=effectiveHabitats(i);
			
			idx_females=find(and(population(habitat,:)==localHabitat, population(sex,:)==0));
			femalesOnHabitatsCell{i} = idx_females;
			femaleFitness=FitnessLite(mean(population(envTraitLocus,idx_females)),habitatCenters(localHabitat),fitnessIOS);
			femalesFitnessOnHabitatsCell{i}=femaleFitness;
			
			idx_males=find(and(population(habitat,:)==localHabitat, population(sex,:)==1));
			malesOnHabitatsCell{i} = idx_males;
			maleFitness=FitnessLite(mean(population(envTraitLocus,idx_males)),habitatCenters(localHabitat),fitnessIOS);
			malesFitnessOnHabitatsCell{i}=maleFitness;
		end
		
		meanFemaleFitness = mean(cell2mat(femalesFitnessOnHabitatsCell));
		meanMaleFitness = mean(cell2mat(malesFitnessOnHabitatsCell));
		
		asexual_females_fitness = Fitness(population(habitat,asexual_females), mean(population(envTraitLocus,asexual_females)), habitatCenters, fitnessIOS);
		if length(asexual_females_fitness)==0
			meanAsexualFemaleFitness=NaN;
		else
			meanAsexualFemaleFitness = mean(asexual_females_fitness);
		end
		%% 6,7,8 ->Mean dispersal probability for sexual females, males, and asexual females: meanProbFemaleDispersal, meanProbMaleDispersal, meanProbAsexualFemaleDispersal
		idx_females = cell2mat(femalesOnHabitatsCell);
		femaleDispersalTraitVector = mean(population(dispLociF(1):dispLociS(end),:))(idx_females);
		femaleDispersalProbability = 1./(1+exp(-femaleDispersalTraitVector));
		meanProbFemaleDispersal = mean(femaleDispersalProbability);
		
		idx_males = cell2mat(malesOnHabitatsCell);
		maleDispersalTraitVector = mean(population(dispLociS(1):dispLociM(end),:))(idx_males);
		maleDispersalProbability = 1./(1+exp(-maleDispersalTraitVector));
		meanProbMaleDispersal = mean(maleDispersalProbability);
		
		if length(asexual_females_fitness)==0
			meanProbAsexualFemaleDispersal=NaN;
		else
			asexualFemaleDispersalTraitVector = mean(population(dispLociF(1):dispLociS(end),:))(asexual_females);
			asexualFemaleDispersalProbability = 1./(1+exp(-asexualFemaleDispersalTraitVector));
			meanProbAsexualFemaleDispersal = mean(asexualFemaleDispersalProbability);
		end
		
		%% Record in SummaryData
		SummaryData(:,t)=[totalNumIndividuals; totalEffectiveHabitats; meanFemaleFitness; meanMaleFitness; meanAsexualFemaleFitness; meanProbFemaleDispersal; meanProbMaleDispersal; meanProbAsexualFemaleDispersal; freqAsexualFemales];
		%%%%%%%%%%%%%%%%%%%%%  end population census  %%%%%%%%%%%%%%%%%%%%%%%
		
		
		%%%%%%%%%%%%%%%%%%%%% start producing next generation %%%%%%%%%%%%%%%%%%%%%
		young=NaN(numLoci,popSize); 
		
		%% Consider a world where selection on Females are hard, and selection on males are soft.
		%% Females compete with each other globally, while males compete with other males on the same habitat.
		for i=1:popSize 
			% pick a mom... 
			sexual_female_fitness=cell2mat(femalesFitnessOnHabitatsCell);
			chosenMom = myID_parent([sexual_female_fitness, asexual_females_fitness]);
			MomID = [idx_females,asexual_females](chosenMom);
			momGenome=population(:,MomID);
			if ismember(MomID,asexual_females)
				young(:,i) = momGenome;
			else
				% pick a dad
				momHabitat = population(habitat,MomID);
				chosenDad=myID_parent(malesFitnessOnHabitatsCell{find(effectiveHabitats==momHabitat)});
				DadID = malesOnHabitatsCell{find(effectiveHabitats==momHabitat)}(chosenDad);
				dadGenome=population(:,DadID);
				%construct offspring genome, assume all loci are independent in linkage
				young(1,i)=momGenome(unidrnd(2)); %Environmental trait from mom
				young(2,i)=dadGenome(unidrnd(2)); %Environmental trait from dad
				young(3,i)=momGenome(unidrnd(2) + 2); %Females dispersal trait from mom
				young(4,i)=dadGenome(unidrnd(2) + 2); %Females dispersal trait from dad
				young(5,i)=momGenome(unidrnd(2) + 4); %shared dispersal trait from mom
				young(6,i)=dadGenome(unidrnd(2) + 4); %shared dispersal trait from dad
				young(7,i)=momGenome(unidrnd(2) + 6); %male dispersal trait from mom
				young(8,i)=dadGenome(unidrnd(2) + 6); %male dispersal trait from dad
				young(sex,i)=unidrnd(2)-1; %assume balanced sex ratio
				young(habitat,i)=momGenome(habitat);
				young(probMating,i)=momGenome(probMating);
			end
			
			%offspring mutation
			mutations=[normrnd(0,sigmaTraitCoordinates,2,1); normrnd(0,sigmaDispersal,6,1); 0; 0; normrnd(0,sigmaTraitCoordinates,1,1)];
			young(:,i)=young(:,i)+mutations;
			young(probMating,:)=max(min(young(probMating,:),1),0); %making sure the probability is a probability

		end
		
		population=young;
		%%%%%%%%%%%%%%%%%%%%%% end producing next generation %%%%%%%%%%%%%%%%%%%%%%
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start dispersal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		traitsOffspringDispersal = mean(population(dispLociF(1):dispLociS(end),:)).*(population(sex,:)==0)+mean(population(dispLociS(1):dispLociM(end),:)).*(population(sex,:)==1)+mean(population(dispLociF(1):dispLociS(end),:)).*(population(sex,:)==2);
		probOffspringDispersal = 1./(1+exp(-traitsOffspringDispersal));
		population(habitat,:) = Disperse(numHabitats, population(habitat,:), probOffspringDispersal, dispersalDeathRate);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end dispersal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start environment change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		envrionmentVector=UpdateEnvironment(envrionmentVector,numGrid,pSpatial,pTemporal,1);
		habitatCenters=envrionmentVector(habitatLocations);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end environment change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
	
	%%%%%%%%%%%%%%%%%%%% Write SummaryData to Output %%%%%%%%%%%%%%%%%%%%%%%%%%%
	csvwrite([strcat("ToMateOrNot-Monogamy-InitDispProb-", num2str(init_Prob_dispersal), "-numGrid-", num2str(numGrid), "-numHabitats-", num2str(numHabitats), "-Rep-", num2str(numRepeats), ".csv")], SummaryData);
end