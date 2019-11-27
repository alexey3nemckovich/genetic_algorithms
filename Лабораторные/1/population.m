classdef population
    
    properties
        % Population properties
        Size = 0
        ChromosomeCount = 0
        GenCount = 0
        MaxGenerationNumber = 0
        MinTValue = -10
        MaxTValue = 10
        
        % Population state
        GenerationNumber = 0
        Population = []
        BestIndivid
        MedianFValue
    end
    
    methods
        function obj = population(size, chromosomeCount, GenCount, maxGenerationNumber, minTValue, maxTValue)
            % set properties
            obj.Size = size;
            obj.ChromosomeCount = chromosomeCount;
            obj.GenCount = GenCount;
            obj.MaxGenerationNumber = maxGenerationNumber;
            obj.MinTValue = minTValue;
            obj.MaxTValue = maxTValue;
            
            obj = obj.reset();
            
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
            obj.GenerationNumber = 0;
        end
        
        function obj = doSelection(obj)
            % selection by ranging
            
            [~, ind] = sort([obj.Population.FValue], 'descend');
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
            % single point crossingover
            
            i = 1;
            while (i+1) <= obj.Size
                if(rand() <= p)
                    r = randi([1 obj.GenCount-1]);
                    r = r + 1; %copying from next
                    
                    for j=1:obj.ChromosomeCount
                        temp = obj.Population(i).Chromosomes(j, r:end);
                        obj.Population(i).Chromosomes(j, r:end) = obj.Population(i+1).Chromosomes(j, r:end);
                        obj.Population(i+1).Chromosomes(j, r:end) = temp;
                    end
                end
                
                i=i+2;
            end
        end
        
        function obj = doMutation(obj, p)
            % binary mutation
            
            for i=1:obj.Size
                for j=1:obj.ChromosomeCount
                    if(rand() <= p)
                        r = randi([1 obj.GenCount]);
                        
                        obj.Population(i).Chromosomes(j, r) = 1 - obj.Population(i).Chromosomes(j, r);
                    end
                end
            end
        end
        
        function obj = calcPopulationParameters(obj)
            bestIndividIndex = 1;
            maxFValue = 0;
            
            for i=1:obj.Size
                
                obj.Population(i).XValues = [];
                
                for j=1:obj.ChromosomeCount
                    doubleXValue = obj.chromosomeToDouble(obj.Population(i).Chromosomes(j,:));
                    obj.Population(i).XValues = [obj.Population(i).XValues, doubleXValue]; 
                end
                
                obj.Population(i).FValue = obj.fitnessValue(obj.Population(i).XValues);
                
                if (i == 1) || (obj.Population(i).FValue > maxFValue)
                    maxFValue = obj.Population(i).FValue;
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
            
            ind.Chromosomes = randi([0,1], obj.ChromosomeCount, obj.GenCount);
        end
        
        function res = chromosomeToDouble(obj, chromosome)
            bin_number_str = num2str(chromosome);
            decimalValue = bin2dec(bin_number_str);
            
            res = obj.MinTValue + decimalValue * (obj.MaxTValue - obj.MinTValue) / (2^obj.GenCount - 1);
        end
        
        function res = fitnessValue(obj, xValues)
            x = xValues(1);
            res = FitnessFunction(x);
        end
        
        function individ = getBestIndivid(obj)
            individ = obj.BestIndivid;
        end
        
    end
end
