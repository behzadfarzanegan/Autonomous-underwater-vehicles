function output = Critic_NL_gamma_bah(x)


%  
%  syms x [1 5]
DS = [];

k=1;
while k <= length(x)
    for i=k:length(x)
      DS =[DS;x(k)*x(i)]; 
    end
    k=k+1;
end
% for k=1:length(x)
%     for j=k:length(x)
%         for l=j:length(x)
%             for i=l:length(x)
%                 DS =[DS;x(k)*x(j)*x(l)*x(i)];
%             end
%         end
%     end
% end

output = DS;


end