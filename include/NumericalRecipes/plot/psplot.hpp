#pragma once
#include <NumericalRecipes/types.hpp>

namespace NumericalRecipes {
namespace psplot{

    namespace {
        using namespace types;
    }

    struct Page {

        static FILE *PLT;
        static char *file;
        char fontname[128];
        double fontsize;

        Page(char *filnam) {
            assign_file(filnam);
        }

        Page(std::string filnam){
            char *cstr = new char[filnam.size()];
            strcpy(cstr, filnam.c_str());
            assign_file(cstr);
        }

        Page() {}

        ~Page() {if (PLT) close();}

        void assign_file(char * filnam){
            file = new char[128];
            strcpy(file,filnam);
            PLT = fopen(file,"wb");
            if (!PLT) throw("failure opening output file for plot");
            fprintf(PLT,"%%!\n/mt{moveto}def /lt{lineto}def /np{newpath}def\n");
            fprintf(PLT,"/st{stroke}def /cp{closepath}def /fi{fill}def\n");
            fprintf(PLT,"/zp {gsave /ZapfDingbats findfont exch ");
            fprintf(PLT,"scalefont setfont moveto show grestore} def\n");
            setfont("Times-Roman",12.);
            setlinewidth(0.5);

        }

        void setfont(char *fontnam, double size) {
            strcpy(fontname,fontnam);
            fontsize = size;
            fprintf(PLT,"/%s findfont %g scalefont setfont\n",fontnam,size);
        }

        void setcolor(int r, int g, int b) {
            fprintf(PLT,"%g %g %g setrgbcolor\n",r/255.,g/255.,b/255.);}

        void setdash(char *patt, int phase=0) {
            fprintf(PLT,"[%s] %d setdash\n",patt,phase);}

        void setlinewidth(double w) {fprintf(PLT,"%g setlinewidth\n",w);}

        void setgray(double w) {fprintf(PLT,"%g setgray\n",w);}

        void gsave() {fprintf(PLT,"gsave\n");}

        void grestore() {fprintf(PLT,"grestore\n");}

        void rawps(char *text) {fprintf(PLT,"%s\n",text);}

        void addtext(char *text) { fprintf(PLT,"(%s) show ",text); }

        void puttext(char *text, double x, double y, double rot=0.0) {
            fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt ",x,y,rot);
            addtext(text);
            fprintf(PLT,"grestore \n");
        }

        void putctext(char *text, double x, double y, double rot=0.0) {
            fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt (%s) ",x,y,rot,text);
            fprintf(PLT,"dup stringwidth pop 2 div neg 0 rmoveto show grestore\n");
        }

        void putrtext(char *text, double x, double y, double rot=0.0) {
            fprintf(PLT,"gsave %g %g translate %g rotate 0 0 mt (%s) ",x,y,rot,text);
            fprintf(PLT,"dup stringwidth pop neg 0 rmoveto show grestore\n");
        }

        void close() {fprintf(PLT,"showpage\n"); fclose(PLT); PLT = NULL;}

        void display() {
            char cmd[128];
            if (PLT) 
                close();
            strcpy(cmd,"\"C:\\Program Files\\gs\\gs9.52\\bin\\gswin64.exe\" ");
            strcat(cmd,file);
            std::cout << cmd << std::endl;
            system(cmd);
        }

        void pointsymbol(double x, double y, int num, double size) {
            fprintf(PLT,"(\\%03o) %g %g %g zp\n",num,x-0.394*size,y-0.343*size,size);
        }

        void lineseg(double xs, double ys, double xf, double yf) {
            fprintf(PLT,"np %g %g mt %g %g lt st\n",xs,ys,xf,yf);
        }

        void polyline(VectorDouble &x, VectorDouble &y, int close=0, int fill=0, int clip=0) {
            int i,n=MIN(x.size(),y.size());
            fprintf(PLT,"np %g %g mt\n",x[0],y[0]);
            for (i=1;i<n;i++) fprintf(PLT,"%g %g lt\n",x[i],y[i]);
            if (close || fill || clip) fprintf(PLT,"cp ");
            if (fill) 
                fprintf(PLT,"fi\n");
            else if (clip) 
                fprintf(PLT,"clip\n");
            else fprintf(PLT,"st\n");
        }
    };

    //---------------------------------------------------------------------------------------------------------

    struct Axis : Page {

        double pll,qll,pur,qur;
        double xll,yll,xur,yur;
        VectorDouble xbox,ybox;
        double majticsz,minticsz;

        Axis(Page &page, double ppll, double ppur, double qqll, double qqur)
            : pll(ppll), qll(qqll), pur(ppur), qur(qqur),
            xll(ppll), yll(qqll), xur(ppur), yur(qqur), xbox(4), ybox(4),
            majticsz(8.), minticsz(4.) {
                strcpy(fontname,page.fontname);
                fontsize = page.fontsize;
                setlimits(xll,xur,yll,yur);
            }

        double p(double x) {return pll + (pur-pll)*(x-xll)/(xur-xll);}
        double q(double y) {return qll + (qur-qll)*(y-yll)/(yur-yll);}
        
        void setlimits(double xxll, double xxur, double yyll, double yyur) {
            xbox[0] = xbox[3] = xll = xxll; ybox[0] = ybox[1] = yll = yyll;
            xbox[1] = xbox[2] = xur = xxur; ybox[2] = ybox[3] = yur = yyur;
        }
        
        void lineseg(double xs, double ys, double xf, double yf) {
            Page::lineseg(p(xs),q(ys),p(xf),q(yf));
        }
        
        void polyline(VectorDouble &x, VectorDouble &y, int close=0, int fill=0, int clip=0) {
            int i;
            VectorDouble xx(x), yy(y);
            for (i=0;i<x.size();i++) xx[i] = p(x[i]);
            for (i=0;i<y.size();i++) yy[i] = q(y[i]);
            Page::polyline(xx,yy,close,fill,clip);
        }
        
        void dot(double x, double y, double size=2.) {
            Page::pointsymbol(p(x),q(y),108,size);
        }
        
        void pointsymbol(double x, double y, int num, double size) {
            Page::pointsymbol(p(x),q(y),num,size);
        }
        
        void lineplot(VectorDouble &x, VectorDouble &y) {polyline(x,y);}
        
        void frame() {polyline(xbox,ybox,1,0);}
        
        void clear() {gsave(); setgray(1.); polyline(xbox,ybox,1,1); grestore();}
        
        void clip() {gsave(); polyline(xbox,ybox,1,0,1);}
        
        void clip(VectorDouble &x, VectorDouble &y) {gsave(); polyline(x,y,1,0,1);}
        
        void unclip() {grestore();}
        
        void xlabel(char *text) {putctext(text,0.5*(pll+pur),qll-2.*fontsize-8.);}
        
        void ylabel(char *text) {putctext(text,pll-3.*fontsize-8.,0.5*(qll+qur),90.);}
        
        void label(char *text, double x, double y, double rot=0.) {puttext(text,p(x),q(y),rot);}
        
        void scalestr(char *str, double x) {
            if (abs(x) < 1.e-15) x = 0.;
            sprintf(str,"%g",x);
        }
        
        void scales(double xmajd, double xmind, double ymajd, double ymind,
        int dox=2, int doy=2, int doxx=1, int doyy=1) {
            char str[128];
            double x,y,xlo,ylo;
            if (dox || doxx) {
                xlo = ceil(MIN(xll,xur)/xmajd)*xmajd;
                for (x=xlo;x<=MAX(xll,xur);x+=xmajd) {
                    scalestr(str,x);
                    if (dox>1) putctext(str,p(x),qll-fontsize-2.);
                    if (dox) Page::lineseg(p(x),qll,p(x),qll+majticsz);
                    if (doxx) Page::lineseg(p(x),qur,p(x),qur-majticsz);
                }
                xlo = ceil(MIN(xll,xur)/xmind)*xmind;
                for (x=xlo;x<=MAX(xll,xur);x+=xmind) {
                    if (dox) Page::lineseg(p(x),qll,p(x),qll+minticsz);
                    if (doxx) Page::lineseg(p(x),qur,p(x),qur-minticsz);
                }
            }
            if (doy || doyy) {
                ylo = ceil(MIN(yll,yur)/ymajd)*ymajd;
                for (y=ylo;y<=MAX(yll,yur);y+=ymajd) {
                    scalestr(str,y);
                    if (doy>1) putrtext(str,pll-4.,q(y)-0.3*fontsize);
                    if (doy) Page::lineseg(pll,q(y),pll+majticsz,q(y));
                    if (doyy) Page::lineseg(pur,q(y),pur-majticsz,q(y));
                }
                ylo = ceil(MIN(yll,yur)/ymind)*ymind;
                for (y=ylo;y<=MAX(yll,yur);y+=ymind) {
                    if (doy) Page::lineseg(pll,q(y),pll+minticsz,q(y));
                    if (doyy) Page::lineseg(pur,q(y),pur-minticsz,q(y));
                }
            }
        }
        
        void autoscales() {
            double xmajd, xmind, ymajd, ymind;
            xmajd = pow(10.,((int)(log10(abs(xur-xll))-1.1)));
            xmind = xmajd/5.;
            ymajd = pow(10.,((int)(log10(abs(yur-yll))-1.1)));
            ymind = ymajd/5.;
            scales(xmajd,xmind,ymajd,ymind);
        }
    };
    FILE *Page::PLT;
    char *Page::file;
    

}//end psplot
}//end NumericalRecipes