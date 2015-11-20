#include "galaxy.h"
#include "star.h"

galaxy::galaxy()
{

}


void galaxy::add(star newstar)
{
    total_stars += 1;
    all_stars.push_back(newstar);
}

