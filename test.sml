structure Test = struct
fun main() =
    let
	open Bfgs
	fun stop _ = OS.Process.exit(OS.Process.success)
	val {b,fail,fmin,fncount,grcount} =
	    vmmin_default(Array.fromList [~1.0,1.0], rosenbrock, rosenbrock1)
	val _ = print("Rosenbrock test expects (1,1):\n" ^
		      "Array.fromList [" ^ Real.toString(Array.sub(b,0)) ^ "," ^
		      Real.toString(Array.sub(b,0)) ^ "]\n")
    in
	stop()
    end
end
val _ = Test.main()
