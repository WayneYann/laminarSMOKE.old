{
	if (minMaxCounter == minMaxUpdate)
	{
		for (label j=0; j<Y.size(); j++)
    		{
			Info << " * Max " << Y[j].name() << " = " << max(Y[j]).value() << endl;		
		}
		
		minMaxCounter = 0;
		
	}

	minMaxCounter++;
}
