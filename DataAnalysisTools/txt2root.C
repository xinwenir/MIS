void PrintVector(TVector3 v)
{
	cout<<"("<<v(0)<<","<<v(1)<<","<<v(2)<<")"<<endl;
}

void SplitString(TString str, TString delim, vector<TString> &vstr){

  vstr.clear();

  TString v;
  Ssiz_t from = 0;
  while (str.Tokenize(v, from, delim)) {
    vstr.push_back(v);
  }

}

void txt2root(string txt_path, string tree_path, double z_aver1, double z_aver2,  double rotAnglex, double rotAngley, int meas_time)
{
	Float_t xmode1,xmode2,ymode1,ymode2;
	Float_t X1,Y1,Z1,X2,Y2,Z2;
	Float_t xo1,zx1,xo2,zx2,yo1,zy1,yo2,zy2;
	Int_t xNBar1,xNBar2,yNBar1,yNBar2;
	Float_t thetax_trk,thetay_trk,thetax,thetay,theta,phi;

	TTree *mytree=new TTree("t_tree","t_tree");

	mytree->Branch("X1",&X1,"X1/F");
	mytree->Branch("Y1",&Y1,"Y1/F");
	mytree->Branch("Z1",&Z1,"Z1/F");
	mytree->Branch("X2",&X2,"X2/F");
	mytree->Branch("Y2",&Y2,"Y2/F");
	mytree->Branch("Z2",&Z2,"Z1/F");

	mytree->Branch("xo1",&xo1,"xo1/F");
	mytree->Branch("zx1",&zx1,"zx1/F");
	mytree->Branch("xo2",&xo2,"xo1/F");
	mytree->Branch("zx2",&zx2,"zx2/F");
	mytree->Branch("yo1",&yo1,"yo1/F");
	mytree->Branch("zy1",&zy1,"zy1/F");
	mytree->Branch("yo2",&yo2,"yo2/F");
	mytree->Branch("zy2",&zy2,"zy2/F");

	mytree->Branch("xNBar1",&xNBar1,"xNBar1/I");
	mytree->Branch("xNBar2",&xNBar2,"xNBar2/I");
	mytree->Branch("yNBar1",&yNBar1,"yNBar1/I");
	mytree->Branch("yNBar2",&yNBar2,"yNBar2/I");

	mytree->Branch("thetax_trk",&thetax_trk,"thetax_trk/F");
	mytree->Branch("thetay_trk",&thetay_trk,"thetay_trk/F");
	mytree->Branch("thetax",&thetax,"thetax/F");
	mytree->Branch("thetay",&thetay,"thetay/F");

	mytree->Branch("theta",&theta,"theta/F");
	mytree->Branch("phi",&phi,"phi/F");

	TTree * tree_txt = new TTree("t_tree", "t_tree");
	tree_txt->ReadFile(txt_path.c_str());

	TFile *f=new TFile(tree_path.c_str(),"RECREATE");

    Float_t x1,y1,x2,y2;
    tree_txt->SetBranchAddress("x1",&x1);
   	tree_txt->SetBranchAddress("x2",&x2);
   	tree_txt->SetBranchAddress("y1",&y1);
   	tree_txt->SetBranchAddress("y2",&y2);
   	tree_txt->SetBranchAddress("zx1",&zx1);
   	tree_txt->SetBranchAddress("zx2",&zx2);
   	tree_txt->SetBranchAddress("zy1",&zy1);
   	tree_txt->SetBranchAddress("zy2",&zy2);
   	tree_txt->SetBranchAddress("xmode1",&xmode1);
   	tree_txt->SetBranchAddress("xmode2",&xmode2);
   	tree_txt->SetBranchAddress("ymode1",&ymode1);
   	tree_txt->SetBranchAddress("ymode2",&ymode2);
   	
    for (int i = 0; i < tree_txt->GetEntries(); ++i)
    {
    	tree_txt->GetEntry(i);

    	xo1=x1;xo2=x2;yo1=y1;yo2=y2;

		xNBar1 = (int)xmode1;
		xNBar2 = (int)xmode2;
		yNBar1 = (int)ymode1;
		yNBar2 = (int)ymode2;

    	Z1=z_aver1;
    	Z2=z_aver2;
		
    	X1=(Z1-zx1)*(xo1-xo2)/(zx1-zx2)+xo1;
    	X2=(Z2-zx1)*(xo1-xo2)/(zx1-zx2)+xo1;
    	Y1=(Z1-zy1)*(yo1-yo2)/(zy1-zy2)+yo1;
    	Y2=(Z2-zy1)*(yo1-yo2)/(zy1-zy2)+yo1;

    	TVector3 v_det_coord1(X1,Y1,Z1);
    	TVector3 v_det_coord2(X2,Y2,Z2);

    	thetax_trk=atan((xo1-xo2)/(zx1-zx2));
    	thetay_trk=atan((yo1-yo2)/(zy1-zy2));

    	// check X1,X2,Y1,Y2, thetax_trk, thetay_trk ...

    	TVector3 v_det_angle(tan(thetax_trk),tan(thetay_trk),1);
    	//check v_det_angle ...

    	TVector3 v_det_coord=v_det_coord1-v_det_coord2;
    	TVector3 v_car_coord(v_det_coord);

    	//check thetax_trk thetay_trk ...

    	TVector3 v_check = (v_det_angle.Unit()).Cross(v_det_coord.Unit());
    	if(v_check.Mag()>1e-5){
    		cerr<<"v_det_angle != v_det_coord: "<<v_check.Mag()<<endl;
    		PrintVector(v_check);
    		PrintVector(v_det_angle.Unit());
    		PrintVector(v_det_coord.Unit());
    	}


    	if (rotAnglex>1e-10 && rotAngley>1e-10)
    	{
    		cerr<<"rotAnglex>1e-10 && rotAngley>1e-10\n";
    		assert(0);
    	}
    	
    	if (rotAnglex<1e-10) 
    		v_car_coord.RotateX(rotAngley);

    	if (rotAngley<1e-10) 
    		v_car_coord.RotateY(rotAnglex);

    	if (rotAnglex<1e-10){
    		thetax=TMath::ATan(v_car_coord(0)/v_car_coord(2));
    		//double thetay_tmath=TMath::ATan(v_car_coord(1)/v_car_coord(2));
    		thetay=thetay_trk+rotAngley;
    	}
    	if (rotAngley<1e-10){
    		thetay=TMath::ATan(v_car_coord(1)/v_car_coord(2));
    		//double thetax_tmath=TMath::ATan(v_car_coord(0)/v_car_coord(2));
    		thetax=thetax_trk+rotAnglex;
    	}

    	//check theta phi ...
    	theta=v_car_coord.Theta();
    	//TMath::ATan(sqrt(pow(tan(thetax),2.0)+pow(tan(thetay),2.0)));
    	phi=v_car_coord.Phi();    	


    	mytree->Fill();
    }

	mytree->Write();

	TVectorD v_meas_time(1);
	v_meas_time[0] = meas_time;

	v_meas_time.Write("meas_time");

	f->Save();
	f->Close();

}