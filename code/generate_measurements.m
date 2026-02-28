diary data/cmdWindowOutput  % write command line output to file

sizes = zeros(1, 64);
times = zeros(3, 64);
errors = zeros(3, 64);

condnos = zeros(1, 64);

for i = 1:64
  i  % print current index

  % load reference solution
  load(fullfile("data","A_" + num2str(i) + ".mat"), "ref_sol")
  
  % keep track of size
  sizes(i) = size(ref_sol, 1);
  
  % load matrix, with same size as the reference solution
  A = logm_testmats(i, sizes(i));
 
  try
      % calculate and save logm cond no
        [X, L, cond] = logm_frechet(A, eye(size(A)));
  catch
      warning("some error in i=" + num2str(i));
      continue
  end
  

  condnos(i) = cond;
  
  %% calculate logm with matlab builtin function
  
  alltimes = zeros(1, 10);
  allerrors = zeros(1, 10);

  for run = 1:11
      tic
      logm_result = logm(A);
      toc
      if run == 1
          continue
      end
      alltimes(run-1) = toc;
      allerrors(run-1) = forwarderror(logm_result, ref_sol);
  end
    
  % save elapsed time
  times(1, i) = mean(alltimes, "all");
  
  % calculate and save forward error
  errors(1, i) = mean(allerrors, "all");
    
  %% calculate logm with complex diss
  
  alltimes = zeros(1, 10);
  allerrors = zeros(1, 10);

  for run = 1:11
      tic
      diss_result = diss(A, 'complex');
      toc
      if run == 1
          continue
      end
      alltimes(run-1) = toc;
      allerrors(run-1) = forwarderror(diss_result, ref_sol);
  end
    
  % save elapsed time
  times(2, i) = mean(alltimes, "all");
  
  % calculate and save forward error
  errors(2, i) = mean(allerrors, "all");
  

  %% if real matrix, also calculate logm with real diss
  if isreal(A)

      alltimes = zeros(1, 10);
      allerrors = zeros(1, 10);
    
      for run = 1:11
          tic
          diss_result_real = diss(A, 'real');
          toc
          if run == 1
              continue
          end
          alltimes(run-1) = toc;
          allerrors(run-1) = forwarderror(diss_result_real, ref_sol);
      end
        
      % save elapsed time
      times(3, i) = mean(alltimes, "all");
      
      % calculate and save forward error
      errors(3, i) = mean(allerrors, "all");
  end
end

save(fullfile("data", "measurements"), "errors", "sizes", "times", "condnos");
diary off