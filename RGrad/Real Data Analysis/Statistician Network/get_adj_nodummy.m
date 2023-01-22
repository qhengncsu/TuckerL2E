load('IDX_to_authors.mat');
gID=idx;
apM=authorPaperBiadj;

Z=zeros(3607,1);
Z(gID)=1;

pZ=Z'*apM;
id2=find(pZ==2);nid2=length(id2);
id3=find(pZ==3);nid3=length(id3);
id4=find(pZ==4);nid4=length(id4);

n=2*nid2+nid3+4*nid4;

Ast=zeros(n,3);

for i=1:nid2
    tid=id2(i);
    tz=apM(:,tid).*Z;
    ids=find(tz>0);
    Ast(2*i-1,1:2)=ids;Ast(2*i-1,3)=ids(1);
    Ast(2*i,1:2)=ids;Ast(2*i,3)=ids(2);
end

for i=1:nid3
    tid=id3(i);
    tz=apM(:,tid).*Z;
    ids=find(tz>0);
    Ast(2*nid2+i,:)=ids;
end

for i=1:nid4
    tid=id4(i);
    tz=apM(:,tid).*Z;
    ids=find(tz>0);
    Ast(2*nid2+nid3+4*(i-1)+1,:)=ids([1,2,3]);
    Ast(2*nid2+nid3+4*(i-1)+2,:)=ids([1,2,4]);
    Ast(2*nid2+nid3+4*(i-1)+3,:)=ids([2,3,4]);
    Ast(2*nid2+nid3+4*(i-1)+4,:)=ids([1,3,4]);
end

Ast=unique(Ast,'rows');

new_gID=gID;
nid=length(new_gID);
A=zeros(nid,nid,nid);

[nedge,tmp]=size(Ast);
for i=1:nedge
    tid1=find(new_gID==Ast(i,1));
    tid2=find(new_gID==Ast(i,2));
    tid3=find(new_gID==Ast(i,3));
    A(tid1,tid2,tid3)=1;A(tid1,tid3,tid2)=1;
    A(tid2,tid1,tid3)=1;A(tid2,tid3,tid1)=1;
    A(tid3,tid1,tid2)=1;A(tid3,tid2,tid1)=1;
end
save('Adj_nodummy.mat','A');





