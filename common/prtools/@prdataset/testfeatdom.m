%TESTFEATDOM Test feature domains
%
%	   [N,I,J] = TESTFEATDOM(A,K,M)
%
% The feature domains of the dataset A are tested. When given,
% the vector K should contain the indices of the features to
% be tested. Optionally, a vector M can be given which indicates which
% of the objects should be taken for testing.
%
% N = 0 if the dataset values are within the domains to be tested.
% N = 1 if somewhere a value is outside a domain. I and J return
%       the indices to the erroneous objects and features.
% 
%	TESTFEATDOM(A,K,M)
%
% Performs the test and prints an error message if a feature value
% outside a domain is encountered.
