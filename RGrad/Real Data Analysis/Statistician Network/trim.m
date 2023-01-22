function [newU]=trim(U,delta)
newU=U;
newU(newU>delta)=delta;
newU(newU<-delta)=-delta;
end