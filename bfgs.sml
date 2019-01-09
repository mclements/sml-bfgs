structure Bfgs : sig
	      type vmminT = {fail: bool, fmin: real, fncount: int, grcount: int, b: real array}
	      val vmmin : real array * (real array -> real) * (real array -> real array) * int * bool * real * real  -> vmminT
	      val vmmin_default : real array * (real array -> real) * (real array -> real array) -> vmminT
	      val rosenbrock : real array -> real
	      val rosenbrock1 : real array -> real array
	      val fdgradient : (real array -> real) -> real -> real array -> real array
	      val fdhessian : (real array -> real array) -> real -> real array -> real Array2.array
	  end = struct

structure Arr = struct
infix 9 sub
open Array
fun dot (v, u) = (
    if length(v) <> length(u) then
	raise ListPair.UnequalLengths
    else ();
    foldli (fn (i, v_i, sum) => Real.+(Real.*(v_i, u sub i), sum)) 0.0 v
)
fun kv(k,v) = tabulate (length v, fn i => Real.*(k,v sub i))
fun vk(v,k) = tabulate (length v, fn i => Real.*(k,v sub i))
fun map2 f v1 v2 =
    tabulate (length v1, fn i => f(v1 sub i, v2 sub i))
fun map f v = tabulate (length v, fn i => f(v sub i))
fun v1+v2 = map2 Real.+ v1 v2
fun v1-v2 = map2 Real.- v1 v2
fun hadamard(v1,v2) = map2 Real.* v1 v2
fun ~v = map Real.~ v
fun print a = 
    let
	val n = length a
	fun loop(i,x) =
	    (if i=0 then TextIO.print "Array.fromList([" else ();
	     if i>0 then TextIO.print "," else ();
	     (TextIO.print o Real.toString) x;
	     if i=Int.-(n,1) then TextIO.print "])\n" else ())
    in
	appi loop a
    end
end

structure Vec = struct
infix 9 sub
open Vector
fun dot (v, u) = (
    if length(v) <> length(u) then
	raise ListPair.UnequalLengths
    else ();
    foldli (fn (i, v_i, sum) => Real.+(Real.*(v_i, u sub i), sum)) 0.0 v
)
fun kv(k,v) = tabulate (length v, fn i => Real.*(k,v sub i))
fun vk(v,k) = tabulate (length v, fn i => Real.*(k,v sub i))
fun map2 f v1 v2 =
    tabulate (length v1, fn i => f(v1 sub i, v2 sub i))
fun map f v = tabulate (length v, fn i => f(v sub i))
fun v1+v2 = map2 Real.+ v1 v2
fun v1-v2 = map2 Real.- v1 v2
fun hadamard(v1,v2) = map2 Real.* v1 v2
fun ~v = map Real.~ v
end

structure Mat = struct
open Array2
fun identity n =
    tabulate RowMajor (n, n, fn (j,k) => if j=k then 1.0 else 0.0)
fun mk(m: real array,k: real) =
    tabulate RowMajor (nRows m, nCols m, fn (i,j) => Real.*(k,sub(m,i,j)))
fun km(k: real, m: real array) =
    tabulate RowMajor (nRows m, nCols m, fn (i,j) => Real.*(k,sub(m,i,j)))
fun mv(m: real array, v: real Array.array) =
    Array.tabulate(nRows m, fn j => Vec.dot(row(m,j),Arr.vector(v)))
fun (m1: real array) * (m2: real array) =
    tabulate RowMajor (nRows m1, nCols m2, fn (i,j) => Vec.dot(row(m1,i),column(m2,j)))
fun map2 f (m1: real array) (m2: real array) =
    tabulate RowMajor (nRows m1, nCols m1, fn (i,j) => f(sub(m1,i,j), sub(m2,i,j)))
fun map f (m: real array) =
    tabulate RowMajor (nRows m, nCols m, fn (i,j) => f(sub(m,i,j)))
fun (m1: real array) + (m2: real array) = map2 Real.+ m1 m2
fun (m1: real array) - (m2: real array) = map2 Real.- m1 m2
fun hadamard(m1: real array,m2: real array) = map2 Real.* m1 m2
fun toList(m: real array) =
    List.tabulate(nRows m, fn i => List.tabulate(nCols m, fn j => sub(m,i,j)))
fun ~m = map Real.~ m
fun print a =
    let
	val (m,n) = dimensions a
	fun loop(i,j,x) =
	    (if i=0 andalso j=0 then TextIO.print "Array2.fromList([" else ();
	     if j=0 andalso i>0 then TextIO.print "," else ();
	     if j=0 then TextIO.print "[" else ();
	     if j>0 then TextIO.print "," else ();
	     (TextIO.print o Real.toString) x;
	     if j=Int.-(n,1) then TextIO.print "]" else ();
	     if i=Int.-(m,1) andalso j=Int.-(n,1) then TextIO.print "])\n" else ())
    in
	appi RowMajor loop {base=a,nrows=SOME m,ncols=SOME n,row=0,col=0}
    end
end

		    
(* infix 8 sub *)
type vmminT = {fail: bool, fmin: real, fncount: int, grcount: int, b: real array}
fun vmmin(b: real array, fminfn: real array -> real, fmingr: real array -> real array,
	  maxit: int, trace: bool, abstol: real, reltol: real) : vmminT =
    if maxit <= 0 then 
        {fail=false, fmin=fminfn b, fncount = 0, grcount = 0, b = b} : vmminT
    else let
	fun diag (i:int,j:int,x:real) = if i=j then 1.0 else 0.0
	and id x = x
	and counting v = Array.foldl (fn (x,agg) => if x then agg+1 else agg) 0 v
	and mcopy m = Mat.map id m
	and acopy a = Arr.map id a
	and incr (x : int ref) = (x := !x+1)
	and printTrace(label,x) = if trace then print(label ^ ": " ^ Real.toString(x) ^ "\n") else ()
	val n = Arr.length b
        val reltest = 10.0
        and stepredn = 0.2
        and acctol = 0.0001
	and mB = Mat.identity n
	and fmin = ref (fminfn b)
        and g = ref (fmingr b)
        and funcount = ref 1
        and gradcount = ref 1
	and iter = ref 1
	and ilast = ref 1
	and bref = ref b
        fun loop1() =
	    let
		val _ = if !ilast = (!gradcount) then (Mat.modifyi Mat.RowMajor diag {base=mB,row=0,col=0,nrows=SOME n,ncols=SOME n}) else ()
                and x = ref(acopy (!bref))
                and c = ref(acopy (!g))
		and count = ref 0
                val t = ref (Arr.~(Mat.mv(mB,!g)))
                val gradproj = Arr.dot(!t,!g)
		val _ = printTrace("gradproj", gradproj)
		val _ = if (gradproj < 0.0) (* search direction is downhill *)
			then
			    let
				fun loop2(steplength) =
				    let
					val _ = (bref := Arr.+(!x, Arr.kv(steplength,!t)))
					val _ = (count := counting
							      (Arr.tabulate(Arr.length (!bref),
									    fn i => Real.==(reltest + Arr.sub(!x,i),
											    reltest + Arr.sub(!bref,i)))))
					val f = if !count=n then Real.posInf else (incr funcount; fminfn (!bref))
					val accpoint = !count < n andalso
						       f<Real.posInf andalso
						       (f <= !fmin + gradproj * steplength * acctol)
				    in
					if !count=n orelse accpoint then (steplength,f) else loop2(steplength*stepredn)
				    end
				val (steplength,f) = loop2(1.0)
				val enough = f > abstol andalso
					     Real.abs(f-(!fmin))>reltol*(Real.abs(!fmin)+reltol)
				val _ = if not enough then (count := n; fmin := f) else ()
			    in
				if !count<n then
				    let
					val _ = (fmin := f;
						 g := fmingr (!bref);
						 incr gradcount;
						 incr iter;
						 t := Arr.kv(steplength,(!t):real array);
						 c := Arr.-(!g,!c))
					val d1 = Arr.dot(!t,!c)
				    in
					if d1>0.0 then
					    let val _ = (x := Mat.mv(mB,!c))
						val d2 = Arr.dot(!x,!c)
						val d3 = 1.0+d2/d1
						val _ = Mat.modifyi
							    Mat.RowMajor
							    (fn (i,j,mBij) => mBij+(d3*Arr.sub(!t,i)*Arr.sub(!t,j)-Arr.sub(!x,i)*Arr.sub(!t,j)-Arr.sub(!t,i)*Arr.sub(!x,j))/d1)
							    {base=mB,row=0,col=0,nrows=SOME n,ncols=SOME n}
						val _ = (printTrace("d1",d1);
							 printTrace("x0",Arr.sub(!x,0));
							 printTrace("c0",Arr.sub(!c,0));
							 printTrace("t0",Arr.sub(!t,0));
							 printTrace("g0",Arr.sub(!g,0));
							 if trace then Mat.print mB else ();
							 printTrace("d2",d2);
							 printTrace("d3",d3))
					    in
						()
					    end (* d1<0.0 *)
					else (ilast := !gradcount)
				    end
				else (* no progress *)
				    if !ilast < !gradcount then (ilast := !gradcount; count := 0) else ()
			    end
			else (* uphill search *)
			    (count := 0; if !ilast = !gradcount then count := n else ilast := !gradcount)
		val _ = if !iter<maxit andalso !gradcount - !ilast > 2*n then (ilast := !gradcount) else ()
	    in
		if !iter >= maxit then {fail=true, fmin = (!fmin), fncount = (!funcount), grcount = (!gradcount), b = !bref} : vmminT
		else
		    if !count=n andalso !ilast = !gradcount then {fail=false, fmin = (!fmin), fncount = (!funcount), grcount = (!gradcount), b = !bref} : vmminT
		    else loop1()
	    end
    in
	loop1()
    end
fun vmmin_default(b, fminfn, fmingr) = vmmin(b, fminfn, fmingr, 100, false, Real.negInf, 1.490116e~08)

(* Rosenbrock function *)
(* Axiom: expr := (a-x)^2+b*(y-x^2)^2; D(expr, x); D(expr, y) *)
local
    val a = 1.0
    val b = 100.0
in
fun rosenbrock par =
    let
	val x = Array.sub(par,0)
	val y = Array.sub(par,1)
    in
	(a-x)*(a-x)+b*(y-x*x)*(y-x*x)
    end
and rosenbrock1 par =
    let
	val x = Array.sub(par,0)
	val y = Array.sub(par,1)
    in
	Array.fromList [~4.0*b*x*y+4.0*b*x*x*x+2.0*x-2.0*a,
		       2.0*b*y-2.0*b*x*x]
    end
end

(* finite differences *)
fun fdgradient f eps x =
    let
	val n = Array.length x
	fun delta(x, i, eps) = Array.tabulate(n, fn j => Array.sub(x,j) + (if i=j then eps else 0.0))
    in
	Array.tabulate(n, fn i => (f(delta(x,i,eps)) - f(delta(x,i,~eps))) / 2.0 / eps)
    end

fun fdhessian grad eps beta =
    let
	val n = Array.length beta
	fun delta(x, i, eps) = Array.tabulate(n, fn (j) => Array.sub(x,j) + (if i=j then eps else 0.0))
	val lol = List.tabulate(n, fn i =>
				      let val upper = grad(delta(beta,i,eps))
					  val lower = grad(delta(beta,i,~eps))
				      in
					  List.tabulate(n, fn i => (Array.sub(upper,i)-Array.sub(lower,i))/2.0/eps)
				      end)
    in
	Array2.fromList lol
    end

end (* struct *)

(* 
(* Tests *)

open Bfgs
vmmin_default(Array.fromList [~1.0,1.0], rosenbrock, rosenbrock1)

vmmin(Array.fromList [~1.0,1.0], rosenbrock, rosenbrock1,
      100, true, Real.negInf, 1.490116e~08)

*)
