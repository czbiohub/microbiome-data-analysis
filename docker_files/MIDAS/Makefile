NAME=midas
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=sunitjain/$(NAME)
ECR=423543210473.dkr.ecr.us-west-2.amazonaws.com/microbiome/$(NAME)

all: build push run

build:
	docker image build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile .
	docker image tag $(NAME):latest $(REGISTRY):latest
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)
	
push:
	docker image push $(REGISTRY):$(TAG)
	docker image push $(REGISTRY):latest

run:
	docker container run --rm $(NAME):latest

release: build push repo-login
	docker image tag $(NAME):latest $(ECR):latest
	docker image tag $(NAME):$(TAG) $(ECR):$(TAG)
	docker image push $(ECR):$(TAG)

## Helpers

# generate script to login to aws docker repo
CMD_REPOLOGIN := "eval $$\( aws ecr"
ifdef AWS_CLI_REGION
CMD_REPOLOGIN += " --region us-west-2"
endif
CMD_REPOLOGIN += " get-login --no-include-email \)"

# login to AWS-ECR
repo-login: ## Auto login to AWS-ECR unsing aws-cli
	@eval $(CMD_REPOLOGIN)