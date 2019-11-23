using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoltApp
{

    public class BoltAlgorithm
    {
       public const double PI = 3.1415926;

        public double Ljx(double nv, double sh, double numx)//x向螺栓距计算(螺栓列数、行距、序号)
        {
            double Lx = 0.0;
            Lx = -0.5 * (nv - 2 * numx + 1) * sh;
            return Lx;
        }
        public double Ljy(double nh, double sv, double numy)//x向螺栓距计算(螺栓列数、行距、序号)
        {
            double Ly = 0.0;
            Ly = 0.5 * (nh - 2 * numy + 1) * sv;
            return Ly;
        }
        public double Ljxx(double nv, double sh, double numxx)//x向螺栓距受压一侧最外排螺栓距离计算(螺栓列数、行距、序号)
        {
            double Lxx = 0.0;
            Lxx = (nv - numxx) * sh;
            return Lxx;
        }

        public double Ljyy(double nh, double sv, double numyy)//y向螺栓距受压一侧最外排螺栓距离计算(螺栓列数、行距、序号)
        {
            double Lyy = 0.0;
            Lyy = (nh - numyy) * sv;
            return Lyy;
        }

        public double Lj(double rx, double ry)//螺栓距计算(x向螺栓距，y向螺栓距)
        {
            double L;
            L = Math.Sqrt(rx * rx + ry * ry);
            return L;
        }

        public double LNx(double sh, double nv)//x向轴力作用点距离受压一侧最外排螺栓的垂直距离
        {
            double Lnx;
            Lnx = 0.5 * (nv - 1) * sh;
            return Lnx;
        }

        public double LNy(double sv, double nh)//y向轴力作用点距离受压一侧最外排螺栓的垂直距离
        {
            double Lny;
            Lny = 0.5 * (nh - 1) * sv;
            return Lny;
        }

        public double Lfvb(double Dl)
        {
            double lfvb = 0.0;
            if (Dl == 4.6)
            {
                lfvb = 140;//Mpa
            }
            else
                if (Dl == 4.8)
            {
                lfvb = 140;//Mpa
            }
            else
                    if (Dl == 5.6)
            {
                lfvb = 190;//Mpa
            }
            else
                        if (Dl == 8.8)
            {
                lfvb = 320;//Mpa
            }
            return lfvb;
        }

        public double Lfcb(double Dl, double Dm)
        {
            double lfcb = 0.0;
            if (Dl == 4.6 || Dl == 4.8)
            {
                if (Dm == 235)
                {
                    lfcb = 305;//Mpa
                }
                else
                    if (Dm == 345)
                {
                    lfcb = 385;//Mpa
                }
                else
                        if (Dm == 390)
                {
                    lfcb = 400;//Mpa
                }
                else
                            if (Dm == 420)
                {
                    lfcb = 425;//Mpa
                }
                else
                                if (Dm == 460)
                {
                    lfcb = 450;//Mpa
                }
            }
            else
                if (Dl == 5.6 || Dl == 8.8)
            {
                if (Dm == 235)
                {
                    lfcb = 405;//Mpa
                }
                else
                    if (Dm == 345)
                {
                    lfcb = 510;//Mpa
                }
                else
                        if (Dm == 390)
                {
                    lfcb = 530;//Mpa
                }
                else
                            if (Dm == 420)
                {
                    lfcb = 560;//Mpa
                }
                else
                                if (Dm == 460)
                {
                    lfcb = 595;//Mpa
                }
            }
            return lfcb;
        }

        public  double Lftb(double Dl)
        {
            double lftb = 0.0;
            if (Dl == 4.6)
            {
                lftb = 170;//Mpa
            }
            else
                if (Dl == 4.8)
            {
                lftb = 170;//Mpa
            }
            else
                    if (Dl == 5.6)
            {
                lftb = 210;//Mpa
            }
            else
                        if (Dl == 8.8)
            {
                lftb = 400;//Mpa
            }
            return lftb;
        }

        public double P(double d)
        {
            double p;
            if (d <= 12)
                p = 1.25 + (d - 8) * 0.125;
            else
                p = 2 + (d - 16) * 0.125;
            return p;

        }

        public double De(double d, double p)
        {
            double de;
            de = (d - 13 * Math.Sqrt(3) * p / 24);
            return de;
        }

        public double Fc(double Dc)
        {
            double fc = 0.0;
            if (Dc == 15)
            {
                fc = 7.2;
            }
            else
                if (Dc == 20)
            {
                fc = 9.6;
            }
            else
                    if (Dc == 25)
            {
                fc = 11.9;
            }
            else
                        if (Dc == 30)
            {
                fc = 14.3;
            }
            else
                            if (Dc == 35)
            {
                fc = 16.7;
            }
            else
                                if (Dc == 40)
            {
                fc = 19.1;
            }
            else
                                    if (Dc == 45)
            {
                fc = 21.1;
            }
            else
                                        if (Dc == 50)
            {
                fc = 23.1;
            }
            else
                                            if (Dc == 55)
            {
                fc = 25.3;
            }
            else
                                                if (Dc == 60)
            {
                fc = 27.5;
            }
            else
                                                    if (Dc == 65)
            {
                fc = 29.7;
            }
            else
                                                        if (Dc == 70)
            {
                fc = 31.8;
            }
            else
                                                            if (Dc == 75)
            {
                fc = 33.8;
            }
            else
                                                                if (Dc == 80)
            {
                fc = 35.9;
            }
            return fc;
        }

    }
}
