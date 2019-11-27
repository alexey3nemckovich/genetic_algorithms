classdef population
    
    properties
        % Population properties
        Size = 0
        ChromosomeCount = 0
        GenCount = 0
        PossibleGens = []
        Values = []
        
        % Population state
        GenerationNumber = 0
        Population = []
        BestIndivid
        MedianFValue
    end
    
    methods
        function obj = population(size, values)
            obj = obj.reset();
            
            % set properties
            obj.ChromosomeCount = 1;
            obj.Size = size;
            obj.GenCount = length(values);
            obj.PossibleGens = values(:,1)';
            obj.Values = values;
            
            % init population
            obj = obj.initPopulation();
        end
        
        function obj = initPopulation(obj)
            obj.GenerationNumber = 1;
            
            for i=1:obj.Size
                obj.Population = [obj.Population, obj.generateNewIndivid()];
            end
            
            obj = obj.calcPopulationParameters();
        end
        
        function obj = goToNextGeneration(obj, crossingoverProbability, mutationProbability)
            % choosing parents (selection)
            obj = obj.doSelection();
            
            % crossing over
            obj = obj.doCrossingOver(crossingoverProbability);
            
            % mutation
            obj = obj.doMutation(mutationProbability);
            
            %update parameters
            obj = obj.calcPopulationParameters();
            
            obj.GenerationNumber = obj.GenerationNumber + 1;
        end
        
        function obj = reset(obj)
            obj.Population = [];
            obj.PossibleGens = [];
            obj.Values = [];
            obj.GenerationNumber = 0;
        end
        
        function obj = doSelection(obj)
            % selection by ranging
            
            [~, ind] = sort([obj.Population.FValue], 'ascend');
            obj.Population = obj.Population(ind);
            
            a = 1+rand();
            b = 2-a;
            
            for i=1:obj.Size
                obj.Population(i).P = (1/obj.Size)*(a - (a-b)*(i-1)/(obj.Size-1));
            end
            
            selectedPopulation = [];
            
            for i=1:obj.Size
                r = rand();
                
                sum = 0;
                for j=1:obj.Size
                    sum = sum + obj.Population(j).P;
                    
                    if sum >= r
                        selectedPopulation = [selectedPopulation, obj.Population(j)];
                        break;
                    end
                end
            end
            
            obj.Population = selectedPopulation;
        end
        
        function obj = doCrossingOver(obj, p)
            % ordered crossingover
            m = 2;
            
            i = 1;
            while (i+1) <= obj.Size
                for j=1:obj.ChromosomeCount
                    if(rand() <= p)
                        r = [];
                        for k=1:m
                            temp = randi([1 obj.GenCount]);
                            while ismember(temp, r) 
                                temp = randi([2 obj.GenCount]);
                            end

                            r = [r, temp];
                        end
                        
                        r = sort(r);
                        
                        leftBorder = r(1);
                        rightBorder = r(2)-1;
                        
                        p1 = obj.Population(i).Chromosomes(j,:);
                        p2 = obj.Population(i+1).Chromosomes(j,:);
                        
                        p1unique = p1(j, leftBorder:rightBorder);
                        p2unique = p2(j, leftBorder:rightBorder);
                        
                        diff1 = p2(rightBorder+1:obj.GenCount);
                        diff1 = [diff1, p2(1:rightBorder)];
                        
                        diff2 = p1(rightBorder+1:obj.GenCount);
                        diff2 = [diff2, p1(1:rightBorder)];
                        
                        diff1 = setdiff(diff1, p1unique,'stable');
                        diff2 = setdiff(diff2, p2unique,'stable');
                        diffLength = length(diff1);
                        
                        rightPartLength = obj.GenCount - rightBorder;
                        leftPartLength = leftBorder - 1;
                        
                        obj.Population(i).Chromosomes(j, rightBorder+1:obj.GenCount) = diff1(1:rightPartLength);
                        obj.Population(i).Chromosomes(j, 1:(leftBorder-1)) = diff1(rightPartLength+1:diffLength);
                        obj.Population(i+1).Chromosomes(j, rightBorder+1:obj.GenCount) = diff2(1:rightPartLength);
                        obj.Population(i+1).Chromosomes(j, 1:(leftBorder-1)) = diff2(rightPartLength+1:diffLength);
                    end
                end
                
                i=i+2;
            end
        end
        
        function obj = doMutation(obj, p)
            % inversion mutation
            
            for i=1:obj.Size
                for j=1:obj.ChromosomeCount
                    if(rand() <= p)
                        r = randi([1 obj.GenCount]);
                        
                        r1 = r;
                        while r1 == r
                            r1 = randi([1 obj.GenCount]);
                        end
                        
                        if r > r1
                            ttemp = r;
                            r = r1;
                            r1 = ttemp;
                        end
                        
                        temp = obj.Population(i).Chromosomes(j, r:r1);
                        temp = flip(temp);
                        obj.Population(i).Chromosomes(j, r:r1) = temp;
                    end
                end
            end
        end
        
        function obj = calcPopulationParameters(obj)
            bestIndividIndex = 1;
            minFValue = 0;
            
            for i=1:obj.Size
                
                obj.Population(i).FValue = obj.fitnessValue(obj.Population(i).Chromosomes);
                
                if (i == 1) || (obj.Population(i).FValue < minFValue)
                    minFValue = obj.Population(i).FValue;
                    bestIndividIndex = i;
                end
            end
            
            if obj.Size > 0
                obj.BestIndivid = obj.Population(bestIndividIndex);
                obj.MedianFValue = median([obj.Population.FValue]);
            end
        end
        
        function ind = generateNewIndivid(obj)
            ind = individ();
            
            ind.Chromosomes = [ind.Chromosomes, obj.PossibleGens(randperm(length(obj.PossibleGens)))];
        end
        
        function res = fitnessValue(obj, chromosome)
            res = 0;
            
            for i=1:obj.GenCount
                
                currTown = obj.Values(chromosome(i),:);
                nextTown = [];
                if i == obj.GenCount
                    nextTown = obj.Values(chromosome(1),:);
                else
                    nextTown = obj.Values(chromosome(i+1),:);
                end
                
                x1 = currTown(2);
                y1 = currTown(3);
                x2 = nextTown(2);
                y2 = nextTown(3);
                
                res = res + sqrt((x1-x2)^2+(y1-y2)^2);
            end
        end
        
        function individ = getBestIndivid(obj)
            individ = obj.BestIndivid;
        end
        
    end
end
