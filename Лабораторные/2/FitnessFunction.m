function res = FitnessFunction(x1, x2)
    res = ((x1.^2+x2.^2).^(1/4)).*(sin(50.*((x1.^2+x2.^2).^(1/10))).^2+1);
end
