% Requires the Advanpix multiprecision computing toolbox.
% the following line (must be adapted) adds it to the path
addpath('C:\Users\Andreas\Documents\Multiprecision Computing Toolbox\')

for i = 1:64
    size_choices = [50 100 150 200];
    % sample choice at random
    size_choice = size_choices(randi([1 4], 1));
    A = logm_testmats(i, size_choice);
    
    % set precision to 64 digits
    mp.Digits(64);
    A_hiprec = mp(A);

    % Perform matrix logarithm on the high precision matrix
    ref_sol = logm(A_hiprec);
    
    ref_sol = double(ref_sol);
    save(fullfile("data", "A_" + num2str(i) + ".mat"), "ref_sol");
end