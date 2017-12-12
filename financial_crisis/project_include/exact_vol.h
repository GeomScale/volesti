




double vol_Ali(std::vector<double> plane, double zit, int dim){
	
	double vol;
	int i,j,J,counter,K,k;
	//double min,max,step,z,zit;
	
	//std::cout<<"hello"<<std::endl;
	//min=min_coeff(plane);
	//max=max_coeff(plane);
	//step=(max-min)/100;
	//std::cout<<"helloooo"<<std::endl;
	double *Y = (double *)malloc((dim+2) * sizeof(double));
	double *X = (double *)malloc((dim+2) * sizeof(double));
	double *a = (double *)malloc((dim+2) * sizeof(double));
	
	//std::cout<<min<<" "<<max<<" "<<step<<" "<<plane.size()<<std::endl;
	
	//z=min;
	//while (z<=max){
		
	J=0; K=0; counter=0;
	if (zit<0){
		X[0]=zit; J++;
	}else{
		Y[0]=zit; counter++;
	}
	
	for (i=0; i<dim; i++){
		
		a[i]=0.0;
		
		if (plane[i]+zit<0){
			X[J]=plane[i]+zit;
			J++;
		}else{
			Y[counter]=plane[i]+zit;
			counter++;
		}
	}
	K=dim+1-J;
	a[0]=1.0; a[dim]=0.0; a[dim+1]=0.0;
	
	for (i=0; i<J; i++){
		for (k=1; k<K+1; k++){
			a[k]=( Y[k-1]*a[k] - X[i]*a[k-1] ) / ( Y[k-1]-X[i] );
		}
	}
		
	//results.push_back(a[K]);
	//z+=step;
	
	vol=a[K];	
	
	free(Y); free(X); free(a);
	return vol;
	
	
}



//Lawrence for computing the volume between two parallel hyperplanes
double cont_lawr_vols(std::vector<double> plane_arx, int dim, double zz1, double zz2){
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// the random engine with this seed
	RNGType rng(seed);
	//RNGType rng;
	//boost::random::uniform_int_distribution<> uidist1(0,dim-1);
	boost::random::uniform_real_distribution<> urdist2(-1.0,1.0);
	int i,j,J,k;
	NT2 max_coef=-1100.0,temp,denom,numer,z1,z2,vol,num,volz2;
	std::vector<NT2> numerators1,numerators2,c,plane,Zs;
	std::vector<double> vols;
	std::vector<int> lefts,lefts_it,rights,rights_it,middles,lefts2,rights2;
	std::pair< std::vector<int>, std::vector<int> > lr,lr2;
	std::pair<NT2,NT2> denoms;
	//c=get_random_hyp(int dim);
	c.clear();
	for (i=0; i<dim; i++){
		plane.push_back(NT2(plane_arx[i]));
		temp=urdist2(rng);
		if (temp>max_coef){
			max_coef=temp;
		}
		c.push_back(NT2(temp));
	}
	//for (i=0; i<Zs2.size(); i++){
	//	Zs.push_back(NT2(Zs2[i]));
	//}
	c.push_back(NT2(0.0));
	//std::cout<<"c is: "<<std::endl;
	//for (i=0; i<dim+1; i++){
	//	std::cout<<c[i]<<std::endl;
	//}
	
	//double q=c[dim],vol,z1,z2,denom,numer;
	z1=NT2(zz1);
	z2=NT2(zz2);
	//double vol2=vol_Ali(plane_arx, -CGAL::to_double(z1), dim);
	//double vol22=vol_Ali(plane_arx, -CGAL::to_double(z2), dim);
	//std::cout<<"vol2 is: "<<vol2<<std::endl;
	//std::cout<<"plane is: "<<std::endl;
	//for (i=0; i<dim; i++){
	//	std::cout<<plane[i]<<std::endl;
	//}
	//std::cout<<"all z are: "<<std::endl;
	//for (i=0; i<Zs.size(); i++){
	//	std::cout<<Zs[i]<<std::endl;
	//}
	//std::cout<<"z is: "<<z1<<std::endl;
	lr=get_left_right_hyp(dim, plane, z1);
	lefts=lr.first; rights=lr.second;
	
	lr2=get_left_right_hyp(dim, plane, z2);  //for second hyp
	lefts2=lr2.first; rights2=lr2.second;
	
	//std::vector<std::vector<int>> edges(dim);
	std::vector<std::vector<NT2>> quant_vert(dim+1);
	
	//std::cout<<"lefts are: "<<std::endl;
	//for (i=0; i<lefts.size(); i++){
		//std::cout<<lefts[i]<<std::endl;
	//}
	//std::cout<<"rights are: "<<std::endl;
	//for (i=0; i<rights.size(); i++){
		//std::cout<<rights[i]<<std::endl;
	//}
	
	for(i=0; i<dim; i++){
		quant_vert[i].resize(dim+1);
	}
	//std::cout<<"hello1"<<std::endl;
	for(i=0; i<dim+1; i++){
		for (j=i+1; j<(dim+1); j++){
			//std::cout<<i<<" "<<j<<std::endl;
			//quant_vert[i][j]=get_denominator(i,j,dim,plane,c,true);   //exei la8os
			quant_vert[i][j]=NT2(-1);
			//std::cout<<"denom for "<<i<<" "<<j<<" is: "<<quant_vert[i][j-1]<<std::endl;
		}
	}
	
	//std::cout<<"hello2"<<std::endl;
	vol=NT2(0.0);
	volz2=NT2(0.0);
	for (i=0; i<lefts.size(); i++){
		for (j=0; j<rights.size(); j++){
			if (rights[j]>lefts[i]){
				//denom=quant_vert[lefts[i]][rights[j]];
				denom=get_denominator(lefts[i],rights[j],dim,plane,c,true);
				quant_vert[lefts[i]][rights[j]]=denom;
				//std::cout<<"denom for "<<lefts[i]<<" "<<rights[j]<<" is: "<<denom<<std::endl;
			}else{
				//denom=quant_vert[rights[j]][lefts[i]];
				denom=get_denominator(lefts[i],rights[j],dim,plane,c,true);
				quant_vert[rights[j]][lefts[i]]=denom;
				//std::cout<<"denom for "<<lefts[i]<<" "<<rights[j]<<" is: "<<denom<<std::endl;
			}
			//denom=get_denominator(lefts[i],rights[j],dim,true);
			//edges[lefts[i]][rights[j]]=denom;
			numer=get_numerator(lefts[i],rights[j],dim,plane,z1,c);
			//std::cout<<"numer for "<<lefts[i]<<" "<<rights[j]<<" is: "<<numer<<std::endl;
			numerators1.push_back(numer);
			//std::cout<<"numer is: "<<numer<<"denom is: "<<denom<<" numer/denom "<<numer/denom<<std::endl;
			vol=vol+numer/denom;
			//std::cout<<"vol is: "<<vol<<std::endl;
		}
		denom=get_denominator(lefts[i],-1,dim,plane,c,false);
		//std::cout<<"denom for "<<lefts[i]<<" is: "<<denom<<std::endl;
		if (lefts[i]==0){
			//num=c[dim];
			numer=NT2(0);
		//	for (k=1; k<dim; k++){
				//numer=std::pow(c[dim],((NT2)dim));
			//	c[dim]=c[dim]*num;
			//}
		}else{
			//numer=std::pow(c[lefts[i]-1]+c[dim],((NT2)dim));
			numer=c[lefts[i]-1];
			num=numer;
			for (k=1; k<dim; k++){
				//numer=std::pow(c[dim],((NT2)dim));
				//c[dim]=c[dim]*num;
				numer=numer*num;
			}
		}
		//std::cout<<"numer is: "<<numer<<"denom is: "<<denom<<" numer/denom "<<numer/denom<<std::endl;
		vol=vol+numer/denom;
	//	std::cout<<"vol is: "<<vol<<std::endl;
	}
	
	//for hyp2
	//std::cout<<"hello3"<<std::endl;
	for (i=0; i<lefts2.size(); i++){
		for (j=0; j<rights2.size(); j++){
			if (rights2[j]>lefts2[i]){
				if (quant_vert[lefts2[i]][rights2[j]]!=NT2(-1)){
					denom=quant_vert[lefts2[i]][rights2[j]];
				}else{
					denom=get_denominator(lefts2[i],rights2[j],dim,plane,c,true);
				}
				//std::cout<<"denom for "<<lefts[i]<<" "<<rights[j]<<" is: "<<denom<<std::endl;
			}else{
				if(quant_vert[rights2[j]][lefts2[i]]!=NT2(-1)){
					denom=quant_vert[rights2[j]][lefts2[i]];
				}else{
					denom=get_denominator(lefts2[i],rights2[j],dim,plane,c,true);
				}
			}
			//std::cout<<"denom for "<<lefts[i]<<" "<<rights[j]<<" is: "<<denom<<std::endl;
			
			//denom=get_denominator(lefts[i],rights[j],dim,true);
			//edges[lefts[i]][rights[j]]=denom;
			numer=get_numerator(lefts2[i],rights2[j],dim,plane,z2,c);
			//std::cout<<"numer for "<<lefts[i]<<" "<<rights[j]<<" is: "<<numer<<std::endl;
			//numerators1.push_back(numer);
			//std::cout<<"numer is: "<<numer<<"denom is: "<<denom<<" numer/denom "<<numer/denom<<std::endl;
			volz2=volz2+numer/denom;
			//std::cout<<"vol is: "<<vol<<std::endl;
		}
		denom=get_denominator(lefts2[i],-1,dim,plane,c,false);
		//std::cout<<"denom for "<<lefts[i]<<" is: "<<denom<<std::endl;
		if (lefts[i]==0){
			//num=c[dim];
			numer=NT2(0);
		//	for (k=1; k<dim; k++){
				//numer=std::pow(c[dim],((NT2)dim));
			//	c[dim]=c[dim]*num;
			//}
		}else{
			//numer=std::pow(c[lefts[i]-1]+c[dim],((NT2)dim));
			numer=c[lefts2[i]-1];
			num=numer;
			for (k=1; k<dim; k++){
				//numer=std::pow(c[dim],((NT2)dim));
				//c[dim]=c[dim]*num;
				numer=numer*num;
			}
		}
		//std::cout<<"numer is: "<<numer<<"denom is: "<<denom<<" numer/denom "<<numer/denom<<std::endl;
		volz2=volz2+numer/denom;
	//	std::cout<<"vol is: "<<vol<<std::endl;
	}
	
	
	
	std::cout<<"vol1 is: "<<CGAL::to_double(vol)<<std::endl;
	std::cout<<"vol2 is: "<<CGAL::to_double(volz2)<<std::endl;
	//std::cout<<"hello3"<<std::endl;
	//std::cout<<"vol2 is: "<<vol2<<std::endl;
	//std::cout<<"vol22 is: "<<vol22<<std::endl;
	
	return (CGAL::to_double(volz2)-CGAL::to_double(vol));
	
}
