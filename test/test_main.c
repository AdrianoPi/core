#include <test.h>

#include <arch/arch.h>

struct stub_arguments {
	int (*test_fnc)(void);
	unsigned tid;
};

static arch_thr_ret_t ARCH_CALL_CONV test_run_stub(void *arg)
{
	struct stub_arguments *args = arg;
	rid = args->tid;
	int ret = args->test_fnc();
	return ret ? ARCH_THR_RET_FAILURE : ARCH_THR_RET_SUCCESS;
}

int main(int argc, char **argv)
{
	(void)argc; (void)argv;
	int ret = 0;
	arch_thr_t threads[test_config.threads_count];
	struct stub_arguments args[test_config.threads_count];

	if(test_config.test_init_fnc && (ret = test_config.test_init_fnc())){
		printf("Test initialization failed with code %d\n", ret);
		return ret;
	}

	for(unsigned i = 0; i < test_config.threads_count; ++i){
		args[i].test_fnc = test_config.test_fnc;
		args[i].tid = i;
		if (arch_thread_create(&threads[i], test_run_stub, &args[i])) {
			return TEST_BAD_FAIL_EXIT_CODE;
		}
	}

	for(unsigned i = 0; i < test_config.threads_count; ++i){
		arch_thr_ret_t thr_ret;
		if (arch_thread_wait(threads[i], &thr_ret)) {
			return TEST_BAD_FAIL_EXIT_CODE;
		}
		if (thr_ret) {
			printf("Thread %u failed the test\n", i);
			return -1;
		}
	}

	if(test_config.test_fini_fnc && (ret = test_config.test_fini_fnc())){
		printf("Test finalization failed with code %d\n", ret);
		return ret;
	}

	return 0;
}