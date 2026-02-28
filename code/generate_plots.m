% load measurements
load(fullfile("data", "measurements.mat"), ...
  "errors", "sizes", "times", "condnos");

% exclude normal matrices
normalmats = zeros(64, 1);
for i = 1:64
    A = logm_testmats(i, 200);
    [U, T] = schur(A, 'complex');  % Generate Schur decomposition
    if isdiag(T)
        normalmats(i) = 1;  % Mark as normal matrix
    end
end

disp(num2str(sum(normalmats)) + "matrices were normal, exluding from " + ...
        "results. Those were the matrices")
disp(find(normalmats))

sizes = sizes(:, ~normalmats);
errors = errors(:, ~normalmats);
times = times(:, ~normalmats);
condnos = condnos(:, ~normalmats);


% print relevant data
sizes


num_measurements = length(sizes);

u = 2^(-53);

% generate error plots here
[condnos_sorted, indices] = sort(condnos, 'descend');

% generate plot for error, sorted by condno

fig1 = figure('Visible','off', 'Position',[0 0 1500 800]);

sc1 = scatter(0:num_measurements-1, errors(1,indices), 500, 'x', 'LineWidth', 500, 'DisplayName','\texttt{logm}','LineWidth', 2);
sc1.LineWidth = 4;
hold on;

sc2 = scatter(0:num_measurements-1, errors(2,indices), 500, 'o', 'LineWidth', 150, 'DisplayName','DISS (complex)','LineWidth', 2);
sc2.LineWidth = 4;

scatter(0:num_measurements-1, errors(3,indices), 300, 'b', 'square', ... 
    'DisplayName', 'DISS (real)','LineWidth', 4);
colororder("reef")
xlabel('$A$ sorted by $\kappa_{\log}(A)u$', 'Interpreter','latex', 'FontSize', 35);
ylabel('Forward error', 'Interpreter', 'latex', 'FontSize',35);
yticks(logspace(-30, 0, 7))
grid on;
set(gca,'yscale','log', 'FontSize', 30)
legend('Interpreter','latex', 'FontSize', 35, 'Location','southeast')

plot(0:num_measurements-1, condnos_sorted*u, 'DisplayName','$\kappa_{\log}(A)u$', 'LineWidth', 2);
ylim([10^(-30) 1])

% Save the first figure
saveas(fig1, fullfile("figures", "errorcomparison.png"));


[~, indices2] = sort(times(1,:), 'descend');

% generate plot for times, sorted by logm time
fig2 = figure('Visible','off', 'Position',[0 0 1500 500]);
scatter(0:num_measurements-1, times(1, indices2), 'x', 'DisplayName','\texttt{logm}', 'LineWidth', 2,'SizeData',90);
hold on;
scatter(0:num_measurements-1, times(2, indices2), 'o', 'DisplayName','DISS (complex)','LineWidth', 2,'SizeData',90);
scatter(0:num_measurements-1, times(3, indices2), 'green', 'square', ...% 'filled', ...
    'DisplayName', 'DISS (real)','LineWidth', 2,'SizeData',90);
xlabel('$A$ sorted by time of \texttt{logm}', 'Interpreter','latex', 'FontSize', 35);
ylabel('Evaluation time (s)', 'Interpreter', 'latex', 'FontSize', 35);
set(gca, 'FontSize', 30)
grid on;
set(gca,'yscale','log');
set(gca, 'xticklabel', [])
legend('Interpreter', 'latex', 'FontSize', 35);
saveas(fig2, fullfile("figures", "compcostabs.png"));


% generate relative plot for times, sorted by logm time
no_real_present = times(3,indices2)==0;
times_diss = times(2, indices2) .* no_real_present + times(3,indices2);  % use real value when available

relative_increase = (times_diss - times(1, indices2))./times(1, indices2);
fig3 = figure('Visible','off', 'Position',[0 0 1500 350]);
clear labels
for i =1:num_measurements
    if relative_increase(i) < 10
        labels(i) = num2str(round(relative_increase(i)));
    else
        labels(i) = "";
    end
end

b = bar(0:num_measurements-1, relative_increase, 'LineWidth', 2, 'FaceColor','flat', Labels=labels, FontSize=20);

for i = 1:num_measurements
    if no_real_present(i)
        b.CData(i,:) = [0.8660 0.3290 0.0000];
    else
        b.CData(i,:) = [0.0660 0.4430 0.7450];
    end
end


xlabel('$A$ sorted by time of \texttt{logm}', 'Interpreter','latex', 'FontSize', 35);
ylabel(sprintf("Rel. increase of\n evaluation time"), 'Interpreter', 'latex', 'FontSize', 35);
set(gca, 'FontSize', 30)
grid on;
set(gca,'yscale','linear');
set(gca, 'xticklabel', []);
yticks([100, 300, 500])

saveas(fig3, fullfile("figures", "compcostrel.png"));

close all