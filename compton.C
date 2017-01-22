
double compton(double x)
{
    TF1 * func = new TF1("","x*(1-(1/(1+2.*x/511)))",0,5000);
    double y = func->Eval(x);
    return y;
}
