#include <iostream>
class MEtObject {

	public:
		MEtObject();

		void setPt(const float pt) { fPt = pt; std::cout << "MEtObject Set(Pt) = " << fPt << std::endl; }
		void setE(const float e) { fE = e; }

		virtual float Pt() const { return fPt; }
		virtual float E() const { return fE;}
		virtual float ECorr(const float fFac) = 0;

	protected:
		float fPt;
		float fE;
};


class JetObject : public MEtObject {
	
	public:
		JetObject();

		float ECorr(const float fFac) const 
		{ 
			std::cout << "JetObject ECorr = " << fE * fFac << std::endl;
			return fE * fFac; 
		}
		

	private:
		int iNtwrs;

};




int main()
{
	MEtObject *m1;
	JetObject j1;
	j1.setPt(5);
	j1.setE(10);
	
	m1 = &j1;

	std::cout << " m1->ECorr() " << m1->ECorr(0.4) << std::endl;

	return 0;
};
