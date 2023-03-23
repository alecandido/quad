///     Purpose:
///         This function maintains the descending ordering in the
///         list of the local error estimated resulting from the
///         interval subdivision process. At each call two error
///         estimates are inserted using the sequential search
///         method, top-down for the largest error estimate and
///         bottom-up for the smallest error estimate.
///
///     Parameters (meaning at output)
///         limit   :   i32
///                     Maximum number of error estimates the list
///                     can contain.
///
///         last    :   i32
///                     Number of error estimates currently in the list.
///
///         maxerr  :   i32
///                     maxerr points to the nrmax-th largest error
///                     estimate currently in the list.
///
///         ermax   :   f64
///                     nrmax-th largest error estimate.
///                     ermax = elist[maxerr-1]
///
///         elist   :   f64
///                     Vector of dimension last containing the error estimates.
///
///         iord    :   Vec<i32>
///                     Vector of dimension last, the first k elements of which
///                     contain pointers to the error estimates, such that
///                     elist(iord(1)),...,  elist(iord(k)) form a decreasing
///                     sequence, with k = last if last <= (limit/2+2), and
///                     k = limit+1-last otherwise.
///
///         nrmax   :   i32
///                     maxerr = iord[nrmax-1].

pub fn qpsrt2(limit : usize, last : usize, maxerr : &mut usize, ermax : &mut f64, elist :& Vec<f64>,iord : &mut Vec<usize>)
             -> (){

    if last <= 2 {
        iord[0] = 1;
        iord.push(2);
        *maxerr = iord[0];
        *ermax = elist[*maxerr-1];
        return()
    }
    //           this part of the routine is only executed if, due to a
    //           difficult integrand, subdivision increased the error
    //           estimate. in the normal case the insert procedure should
    //           start after the nrmax-th largest error estimate.

    let errmax = elist[*maxerr-1];
    let mut i_nrmax = 1;
    let jupbn : usize;
    //  println!("iord pre while: {:?}",iord);

    while i_nrmax > 0 && *ermax > elist[iord[i_nrmax - 1]] && i_nrmax != 1
    {
        iord[i_nrmax] = iord[i_nrmax - 1] ;
        i_nrmax -= 1;
    }

//           compute the number of elements in the list to be maintained
//           in descending order. this number depends on the number of
//           subdivisions still allowed.


    if last > (limit/2+2) { jupbn = limit + 3 - last;}
    else {  jupbn = last;}

    let errmin = elist[last-1];

    //           insert errmax by traversing the list top-down,
    //           starting comparison from the element elist(iord(nrmax+1)).

    let jbnd = jupbn - 1;
    let ibeg = 2;

    //  println!("jupbn:{jupbn},jbnd:{jbnd},ibeg:{ibeg},errmin:{errmin}");

    if ibeg > jbnd {
        iord[jbnd-1] = *maxerr;
        iord.push(last);
        *maxerr = iord[0];
        *ermax = elist[*maxerr-1];
        //  println!("First return");
        return()
    }

    let mut ii : usize = ibeg;
    let mut i_break = 0; // go to work-around
    for i in ibeg..jbnd+1{
        ii = i;
        //  println!("{ii}");
        let isucc = iord[i-1];
        if errmax >= elist[isucc-1] {
            //  println!("{ii} break");
            i_break = 1;
            break
        }
        if i - 1 == last {
            iord.push(isucc);
        }
        else{
            iord[i-2] = isucc;
        }
    }

    if i_break == 0 {
        iord[jbnd-1] = *maxerr;
        iord.push(last);
        *maxerr = iord[0];
        *ermax = elist[*maxerr-1];
        //  println!("Second return");
        return()
    }

    //           insert errmin by traversing the list bottom-up.

    if ii - 1 == last {
        iord.push(*maxerr);
    }
    else{
        iord[ii-2] = *maxerr;
    }
    let mut k = jbnd;

    for _j in ii..jbnd+1{
        let isucc = iord[k-1];
        if errmin < elist[isucc-1] {
            i_break = 1;
            break
        }
        if k == last - 1 {
            iord.push(isucc);
        }
        else{
            iord[k] = isucc;
        }
        k -= 1;
    }

    if i_break == 0 {
        if ii == last{
            iord.push(last);
        }
        else{
            iord[ii-1] = last;
        }
    }
    else{
        if k == last - 1{
            iord.push(last);
        }
        else{
            iord[k] = last;
        }
    }

    *maxerr = iord[0];
    *ermax = elist[*maxerr-1];
}

