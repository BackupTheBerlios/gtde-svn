// combinations --
// Returns the number of combinations of j objects chosen from n objects .
function c = combinations ( n , j )
	c = exp ( gammaln ( n +1) - gammaln ( j +1) - gammaln (n - j +1));
	// If the input where integers , returns also an integer .
	if ( and ( round ( n )== n ) & and ( round ( j )== j ) ) then
		c = round ( c )
	end
endfunction

