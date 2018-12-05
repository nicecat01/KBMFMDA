function everyresult(disease,miRNA,interaction,score)
a=(interaction-1)*(-1);
score=score.*a;
[m,n]=size(score);
for i=1:m
    list=score(i,:);
    listk=list(list~=0);
    [~,k]=size(listk);
    newlist=sort(listk,'descend');
    filename=disease{i};
    str=cell(k,3);
    str(:,1)=disease(i);

    t=find(score(i,:)==newlist(1));
    [~,p]=size(t);
    for x=1:p
        a=miRNA(t(x),1);
        str(x,2)=a;
        str{x,3}=newlist(1);
    end
    for j=2:k
        if newlist(j)==newlist(j-1)
            continue
        end
        t=find(score(i,:)==newlist(j));
        [~,p]=size(t);
        for x=1:p
            a=miRNA(t(x),1);
            str(j-1+x,2)=a;
            str{j-1+x,3}=newlist(j);
        end

    end
%{
    for j=1:k
        t=find(score(i,:)==newlist(j));
        [~,p]=size(t);
        for x=1:p
            a=miRNA(t(x),1);
            str(j-1+x,2)=a;
            str{j-1+x,3}=newlist(j);
        end

    end
%}
    xlswrite(['.\result\everyresult\',filename],str);
end
end