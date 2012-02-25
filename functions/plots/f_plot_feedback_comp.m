function [T Td] = f_plot_feedback_comp(G, C, Cd)
T=feedback(C*G, 1);
Td=feedback(Cd*G, 1);
stepplot(T, Td);
end