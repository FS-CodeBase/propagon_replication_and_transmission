function res = LogL(x,sampling_times,log_Yn,log_cn,cn,data)

LOGL = 0;
for t = 1:length(sampling_times)
    T = sampling_times(t);
    if T == 0
        T = 10^-16; 
    end
    maxGen = ceil(T*60/90);
    gen_vec = cn(0:maxGen,T)/sum(cn(0:maxGen,T));
    gen_vec = (gen_vec>=0.05).*(1:(maxGen+1)); gen_vec(gen_vec==0)=[];
    minGen = min(gen_vec)-1;
    maxGen = max(gen_vec)-1;
    LogCt = rlog(sort(log_cn(minGen:maxGen,T),'descend'));
    aggs_tmp = data{t};
    for d = 1:length(aggs_tmp)
        LOG_YN = [];
        for G = minGen:maxGen
            LOG_YN = [LOG_YN,log_Yn(G,T,aggs_tmp(d),x(1),x(2))];
        end
        LOGL = LOGL + rlog(LOG_YN) - LogCt;
    end
end
res = LOGL;